#include "cito_sqopt.h"

// ***** DESCRIPTION ***********************************************************
// CitoSQOPT class consists of functions that are used to convert the convex
// subproblems used in SCvx to the appropriate form for SQOPT.

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoSQOPT::CitoSQOPT(const mjModel* model) : m(model), cp(model)
{
    // initialize Eigen variables
    desiredPos.resize(6);       deltaPos.resize(6);
    desiredVel.resize(m->nv);   deltaVel.resize(m->nv);
    dKCon.resize(cp.nPair,cp.N);
    // waypoint parameters
    desiredPosWP.resize(6);     deltaPosWP.resize(6);
    tWP = (int) floor(cp.N/2);
    // read task parameters
    YAML::Node params = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/params.yaml");
    std::vector<double> desiredPosInput = { params["desiredFinalPos"].as<std::vector<double>>() };
    std::vector<double> desiredVelInput = { params["desiredFinalVel"].as<std::vector<double>>() };
    std::vector<double> desiredPosWPInput = { params["desiredWayPtPos"].as<std::vector<double>>() };
    desiredPos = Eigen::Map<Eigen::VectorXd>(desiredPosInput.data(), desiredPosInput.size());
    desiredVel = Eigen::Map<Eigen::VectorXd>(desiredVelInput.data(), desiredVelInput.size());
    desiredPosWP = Eigen::Map<Eigen::VectorXd>(desiredPosWPInput.data(), desiredPosWPInput.size());
    controlJointDOF0 = params["controlJointDOF0"].as<int>();
    // read contact model parameters
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    kCon0 = vscm["kCon0"].as<double>();         // upper bound for the virtual stiffness
    // set SQOPT parameters
    nnH     = 6 + m->nv + cp.N*cp.nPair + 6;    // number of non-zero elements of the Hessian
    lencObj = nnH;                              // number of non-zero elements of the linear term
    neA     = cp.N*cp.n*cp.n + (cp.N+1)*cp.n + cp.N*cp.n*cp.m + ((cp.N+1)*cp.n+cp.N*cp.m)*5;
    n       = ((cp.N+1)*cp.n + cp.N*cp.m)*2;    // *2 is for auxiliary variables for l1-norm (trust region)
    nc      = (cp.N+1)*cp.n + ((cp.N+1)*cp.n+cp.N*cp.m)*2 + 1;
    cObj    = new double[lencObj];
    // indices to move
    indMove = new int[nnH];
    // * virtual stiffness variables
    for( int i=0; i<cp.N; i ++ )
    {
        for( int j=0; j<cp.nPair; j++ )
        {
            indMove[i*cp.nPair+j] = (cp.N+1)*cp.n + cp.N*cp.m - 1 - (i*cp.m + j);
        }
    }
    // * final velocity variables for the control joint
    for( int i=0; i<m->nv; i++ )
    {
        indMove[nnH-1-6-6-i] = cp.N*cp.n + m->nv + i;
    }
    // * final pose variables for the control joint
    for( int i=0; i<6; i++ )
    {
        indMove[nnH-1-6-i] = cp.N*cp.n + controlJointDOF0 + i;
    }
    // * waypoint pose variables for the control joint
    for( int i=0; i<6; i++ )
    {
        indMove[nnH-1-i] = tWP*cp.n + controlJointDOF0 + i;
    }
    // initialize & set options for SQOPT
    cvxProb.initialize("", 1);
    cvxProb.setProbName("SubQP");
    cvxProb.setIntParameter("Print level", 0);
    // set the weights
    lenru = 4;      // number of weights
    ru = new double[lenru];
    ru[0] = params["w1"].as<double>();
    ru[1] = params["w2"].as<double>();
    ru[2] = params["w3"].as<double>();
    ru[3] = params["w4"].as<double>();
    cvxProb.setUserR(ru, lenru);
    // set parameters that are dependent on simulation and model
    leniu = 3;      // number of parameters
    iu = new int[leniu];
    iu[0] = cp.nv;
    iu[1] = cp.N;
    iu[2] = cp.nPair;
    cvxProb.setUserI(iu, leniu);
    // setBound parameters
    dUOffset  = (cp.N+1)*cp.n;
    auxOffset = (cp.N+1)*cp.n+cp.N*cp.m;
}
// ***** FUNCTIONS *************************************************************
// qpHx: sets Hx to the H*x part of the quadratic cost to be multiplied by x'
void CitoSQOPT::qpHx(int *nnH, double x[], double Hx[], int *nState,
                     char cu[], int *lencu, int iu[], int *leniu,
                     double ru[], int *lenru)
{
    // get the parameters
    int nv = iu[0], N = iu[1], nPair = iu[2];
    // final x-y position of the control body
    for( int i=0; i<2; i++ )
    {
        // waypoint
        Hx[i] = ru[0]*x[i];
        // final
        Hx[i+6] = ru[0]*x[i+6];
    }
    // final z position and orientation of the control body
    for( int i=2; i<6; i++ )
    {
        // waypoint
        Hx[i] = ru[1]*x[i];
        // final
        Hx[i+6] = ru[1]*x[i+6];
    }
    // final velocity terms
    for( int i=12; i<12+nv; i++ )
    {
        Hx[i] = ru[2]*x[i];
    }
    // stiffness terms
    for( int i=12+nv; i<12+nv+N*nPair; i++ )
    {
        Hx[i] = ru[3]*x[i];
    }
}

// setCObj: sets linear and constant cost terms of the cost
void CitoSQOPT::setCObj(const eigMm X, const eigMd U,
                        double *ru, double *cObj, double &ObjAdd)
{
    // desired change in the final pose and velocity
    deltaPos.setZero(); deltaVel.setZero();
    deltaPos = desiredPos - X.col(cp.N).segment(controlJointDOF0, 6);
    deltaVel = desiredVel - X.col(cp.N).tail(m->nv);
    deltaPosWP.setZero();
    deltaPosWP = desiredPosWP - X.col(tWP).segment(controlJointDOF0, 6);
    // set linear objective terms
    // * x and y position
    for( int i=0; i<2; i++ )
    {
        // waypoint
        cObj[i] = -ru[0]*deltaPosWP(i);
        // final
        cObj[6+i] = -ru[0]*deltaPos(i);
    }
    // * z position and orientation
    for( int i=2; i<6; i++ )
    {
        // waypoint
        cObj[i] = -ru[1]*deltaPosWP(i);
        // final
        cObj[6+i] = -ru[1]*deltaPos(i);
    }
    // * final velocity
    for( int i=0; i<m->nv; i++ )
    {
        cObj[12+i] = -ru[2]*deltaVel(i);
    }
    // 4) virtual stiffness
    double dKConSN = 0;
    for( int i=0; i<cp.N; i++ )
    {
        dKCon.col(i).setZero();
        for( int j=0; j<cp.nPair; j++ )
        {
            dKCon.col(i)[j] = 0-U.col(i)[m->nu+j];
            cObj[12+m->nv+i*cp.nPair+j] = -ru[3]*dKCon.col(i)[j];
        }
        dKConSN += dKCon.col(i).squaredNorm();
    }
    // constant objective term
    ObjAdd = 0.5*(ru[0]*deltaPos.head(2).squaredNorm() +
                  ru[1]*deltaPos.tail(4).squaredNorm() +
                  ru[2]*deltaVel.squaredNorm() +
                  ru[3]*dKConSN);
}

// solveCvx: solves the convex subproblem
void CitoSQOPT::solveCvx(double *xTraj, double r, const eigMm X, const eigMd U,
                         const derTraj Fx, const derTraj Fu, int *isJFree, int *isAFree,
                         double *qposLB, double *qposUB, double *tauLB, double *tauUB)
{
    // fresh start
    auto *indA  = new int[neA];
    auto *locA  = new int[n+1];
    auto *valA  = new double[neA];
    auto *x     = new double[n+nc];
    auto *bl    = new double[n+nc];
    auto *bu    = new double[n+nc];
    auto *pi    = new double[nc];
    auto *rc    = new double[n+nc];
    auto *hs    = new int[n+nc];
    auto *eType = new int[n+nc];
    // initial guess
    for( int i=0; i<n+nc; i++ )
    {
        x[i]  = 0.0;
        bl[i] = -infBnd;  bu[i] = +infBnd;
        hs[i] = 0;  eType[i] = 0;  rc[i] = 0.0;
        if( i>=n ) { pi[i-n] = 0.0; }
    }
    // set linear constraints and bounds
    this->setA(valA, indA, locA, Fx, Fu);
    this->setBounds(r, X, U, bl, bu, isJFree, isAFree, qposLB, qposUB, tauLB, tauUB);
    // sort A and bounds w.r.t. the order in qpHX
    this->sortToMatch(valA, indA, locA, indMove, bl, bu);
    // set linear and constant cost terms
    this->setCObj(X, U, ru, cObj, ObjAdd);
    // solve the convex subproblem
    cvxProb.solve(Cold, this->qpHx, nc, n, neA, lencObj, nnH, iObj,
                  ObjAdd, valA, indA, locA, bl, bu, cObj,
                  eType, hs, x, pi, rc, nS, nInf, sInf, objective);
    // sort x to regular form
    this->sortX(x, indMove);
    // set only the trajectory-related variables (X and U)
    for( int i=0; i<cp.nTraj; i++ )
    {
        xTraj[i] = x[i];
    }
    // cleaning
    delete []indA;  delete []locA; delete []valA;
    delete []x;     delete []bl;   delete []bu;
    delete []pi;    delete []rc;   delete []hs;   delete[]eType;
}

// setBounds: sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
void CitoSQOPT::setBounds(double r, const eigMm X, const eigMd U,
                          double *bl, double *bu, int *isJFree, int *isAFree,
                          double *qposLB, double *qposUB, double *tauLB, double *tauUB)
{
    // decision variables
    // * states
    for( int i=0; i<cp.N+1; i++ )
    {
        for( int j=0; j<m->nv; j++ )
        {
            // ** change in free joint positions: unbounded
            // ** change in joint positions
            if( isJFree[j] == 0 )
            {
                bl[i*cp.n+j] = qposLB[j] - X.col(i)[j];
                bu[i*cp.n+j] = qposUB[j] - X.col(i)[j];
            }
            // ** change in joint velocities: unbounded (already set)
        }
        // * controls
        if( i < cp.N )
        {
            // ** change in joint torques
            for( int j=0; j<m->nu; j++ )
            {
                if( isAFree == 0 )
                {
                    bl[dUOffset+i*cp.m+j] = tauLB[j] - U.col(i)[j];
                    bu[dUOffset+i*cp.m+j] = tauUB[j] - U.col(i)[j];
                }
            }
            // ** change in virtual stiffness
            for( int j=0; j< cp.nPair; j++ )
            {
                bl[dUOffset+i*cp.m+m->nu+j] = 0 - U.col(i)[m->nu+j];
                bu[dUOffset+i*cp.m+m->nu+j] = kCon0 - U.col(i)[m->nu+j];
            }
        }
    }
    // auxiliary variables > 0
    for( int i=auxOffset; i<n; i++ )
    {
        bl[i] = 0;
    }
    // constraints
    // dynamics == 0
    for( int i=n; i<n+(cp.N+1)*cp.n; i++ )
    {
        bl[i] = 0;
        bu[i] = 0;
    }
    // 0 <= absolute value constraints <= inf
    for( int i=n+(cp.N+1)*cp.n; i<n+(cp.N+1)*cp.n+((cp.N+1)*cp.n+cp.N*cp.m)*2; i++ )
    {
        bl[i] = 0;
    }
    // 0 <= trust region <= r
    bl[n+(cp.N+1)*cp.n+((cp.N+1)*cp.n+cp.N*cp.m)*2] = 0;
    bu[n+(cp.N+1)*cp.n+((cp.N+1)*cp.n+cp.N*cp.m)*2] = r;
}

// setA: creates the sparse A matrix for linearized dynamics, auxiliary
// variables, and trust region constraints
void CitoSQOPT::setA(double *valA, int *indA, int *locA, const derTraj Fx, const derTraj Fu)
{
    int colNo = 0, indNo = 0, indTS = 0;

    locA[0] = 0;
    // columns associated with dx[1,...,N]
    for( int i=0; i<cp.N*cp.n; i++ )
    {
        indA[indNo] = i;
        valA[indNo] = -1;
        indNo++;
        // time step index
        indTS = (int) floor(i/cp.n);
        for( int j=0; j<cp.n; j++ )
        {
            indA[indNo] = (indTS+1)*cp.n+j;
            valA[indNo] = Fx.col(indTS)[(j+i*cp.n)%(cp.n*cp.n)];
            indNo++;
        }
        // for auxiliary variables
        indA[indNo] = (cp.N+1)*cp.n + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (cp.N+1)*cp.n + (cp.N+1)*cp.n + cp.N*cp.m + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;
        // column complete
        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with dx[N+1]
    for( int i=0; i<cp.n; i++ )
    {
        indA[indNo] = cp.N*cp.n+i;
        valA[indNo] = -1;
        indNo++;
        // for auxiliary variables
        indA[indNo] = (cp.N+1)*cp.n + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (cp.N+1)*cp.n + (cp.N+1)*cp.n + cp.N*cp.m + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;
        // column complete
        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with du[1,...,N]
    for( int i=0; i<cp.N*cp.m; i++ )
    {
        // time step index
        indTS = (int) floor(i/cp.m);
        for( int j=0; j<cp.n; j++ )
        {
            indA[indNo] = (indTS+1)*cp.n+j;
            valA[indNo] = Fu.col(indTS)[(j+i*cp.n)%(cp.n*cp.m)];
            indNo++;
        }
        // for auxiliary variables
        indA[indNo] = (cp.N+1)*cp.n + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (cp.N+1)*cp.n + (cp.N+1)*cp.n + cp.N*cp.m + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;
        // column complete
        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with auxiliary variables
    int auxNo1 = 0, auxNo2 = 0;
    for( int i=0; i<(cp.N+1)*cp.n+cp.N*cp.m; i++ )
    {
        indA[indNo] = (cp.N+1)*cp.n + auxNo1;
        valA[indNo] = +1.0;
        indNo++; auxNo1++;
        indA[indNo] = (cp.N+1)*cp.n + (cp.N+1)*cp.n + cp.N*cp.m + auxNo2;
        valA[indNo] = +1.0;
        indNo++; auxNo2++;
        // for l1-norm
        indA[indNo] = (cp.N+1)*cp.n + 2*((cp.N+1)*cp.n + cp.N*cp.m);
        valA[indNo] = +1.0;
        indNo++;
        // column complete
        colNo++;
        locA[colNo] = indNo;
    }
}

// sortToMatch: modifies A and the bounds such that non-zero elements in H come first
void CitoSQOPT::sortToMatch(double *valA, int *indA, int *locA, int *moveIndices, double *bl, double *bu)
{
    int less_counter = 0, iMove = 0;
    for( int i=0; i<nnH; i++ )
    {
        if( i>0 && moveIndices[i]<moveIndices[i-1] ) { less_counter++; }
        iMove = moveIndices[i]+less_counter;
        this->moveColA(valA, indA, locA, iMove);
        this->moveRowBounds(bl, bu, iMove);
    }
}

// moveColA: moves column iMove in A to left
void CitoSQOPT::moveColA(double *valA, int *indA, int *locA, int iMove)
{
    int indAtemp[neA], neCol[n]; double valAtemp[neA];

    for( int i=0; i<n; i++ )
    {
        neCol[i] = locA[i+1]-locA[i];
    }
    for( int i=0; i<neCol[iMove]; i++)
    {
        indAtemp[i] = indA[locA[iMove]+i];
        valAtemp[i] = valA[locA[iMove]+i];
    }
    for( int i=0; i<locA[iMove]; i++ )
    {
        indAtemp[neCol[iMove]+i] = indA[i];
        valAtemp[neCol[iMove]+i] = valA[i];
    }

    locA[0] = 0;  locA[1] = neCol[iMove];
    for( int i=1; i<iMove; i++ )
    {
        locA[i+1] = locA[i] + neCol[i-1];
    }
    for( int i=0; i<locA[iMove+1]; i++ )
    {
        indA[i] = indAtemp[i]; valA[i] = valAtemp[i];
    }
}

// moveRowBounds: moves row iMove in bounds to top
void CitoSQOPT::moveRowBounds(double *bl, double *bu, int iMove)
{
    double bLTemp[iMove], bUTemp[iMove];
    bLTemp[0] = bl[iMove];
    bUTemp[0] = bu[iMove];
    for( int i=0; i<iMove+1; i++ )
    {
        if( i<iMove )
        {
            bLTemp[i+1] = bl[i];
            bUTemp[i+1] = bu[i];
        }
        bl[i] = bLTemp[i];
        bu[i] = bUTemp[i];
    }
}

// sortX: sorts decision variables back to original order
void CitoSQOPT::sortX(double *x, int *moveIndices)
{
    int k = 0;
    xTemp = new double[n];
    for( int i=0; i<n; i++ )
    {
        xTemp[i] = x[nnH+k];
        for( int j=0; j<nnH; j++ )
        {
            if( i == moveIndices[j] )
            {
                xTemp[moveIndices[j]] = x[nnH-j-1];
                k = k - 1;
            }
        }
        k = k + 1;
    }
    for( int i=0; i<n; i++ )
    {
        x[i] = xTemp[i];
    }
}