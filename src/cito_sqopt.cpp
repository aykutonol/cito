#include "cito_sqopt.h"

// ***** DESCRIPTION ***********************************************************
// CitoSQOPT class consists of functions that are used to convert the convex
// subproblems used in SCvx to the appropriate form for SQOPT.

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoSQOPT::CitoSQOPT()
{
    // read task parameters
    YAML::Node paramTask = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/task.yaml");
    std::vector<double> desiredPoseInput = { paramTask["desiredFinalPose"].as<std::vector<double>>() };
    std::vector<double> desiredVeloInput = { paramTask["desiredFinalVelo"].as<std::vector<double>>() };
    desiredPose = Eigen::Map<Eigen::Matrix<double,  6, 1>>(desiredPoseInput.data(), desiredPoseInput.size());
    desiredVelo = Eigen::Map<Eigen::Matrix<double, NV, 1>>(desiredVeloInput.data(), desiredVeloInput.size());
    controlJointDOF0 = paramTask["controlJointDOF0"].as<int>();
    // read contact model parameters
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    kCon0 = vscm["kCon0"].as<double>();
    // trajectories
    dKCon.resize(NTS);
    // indices to move
    // *** virtual stiffness variables
    for( int i=0; i < NTS; i ++ )
    {
        for( int j=0; j<NPAIR; j++ )
        {
            indMove[i*NPAIR+j] = (NTS+1)*N + NTS*M - 1 - (i*M + j);
        }
    }
    // *** final velocity variables for the control joint
    for( int i=0; i<NV; i++ )
    {
        indMove[nnH-1-6-i] = NTS*N + NV + i;
    }
    // *** final pose variables for the control joint
    for( int i=0; i<6; i++ )
    {
        indMove[nnH-1-i] = NTS*N + controlJointDOF0 + i;
    }
    // initialize & set options for SQOPT
    cvxProb.initialize("", 1);
    cvxProb.setProbName("SubQP");
    cvxProb.setIntParameter("Print level", 0);
    // set the weights
    ru[0] = paramTask["w1"].as<double>();
    ru[1] = paramTask["w2"].as<double>();
    ru[2] = paramTask["w3"].as<double>();
    ru[3] = paramTask["w4"].as<double>();
    cvxProb.setUserR(ru, lenru);
}
// ***** FUNCTIONS *************************************************************
// qpHx: sets Hx to the H*x part of the quadratic cost to be multiplied by x'
void qpHx(int *nnH, double x[], double Hx[], int *nState,
          char cu[], int *lencu, int iu[], int *leniu,
          double ru[], int *lenru)
{
    // final x-y position of the control body
    for( int i=0; i<2; i++ )
    {
        Hx[i] = ru[0]*x[i];
    }
    // final z position and orientation of the control body
    for( int i=2; i<6; i++ )
    {
        Hx[i] = ru[1]*x[i];
    }
    // final velocity terms
    for( int i=6; i<6+NV; i++ )
    {
        Hx[i] = ru[2]*x[i];
    }
    // stiffness terms
    for( int i=6+NV; i<6+NV+NTS*NPAIR; i++ )
    {
        Hx[i] = ru[3]*x[i];
    }
}

// setCObj: sets linear and constant cost terms of the cost
void CitoSQOPT::setCObj(const stateTraj X, const ctrlTraj U,
                        double *ru, double *cObj, double &ObjAdd)
{
    // set linear objective terms
    // desired change in the final pose and velocity
    deltaPose.setZero(); deltaVelo.setZero();
    for( int i=0; i<6; i++ )
    {
        deltaPose[i] = desiredPose[i] - X[NTS][controlJointDOF0+i];
    }
    for( int i=0; i<NV; i++ )
    {
        deltaVelo[i] = desiredVelo[i] - X[NTS][controlJointDOF0+NV+i];
    }
    // final position
    for( int i=0; i<2; i++ )
    {
        cObj[i] = -ru[0]*deltaPose[i];
    }
    // final orientation
    for( int i=2; i<6; i++ )
    {
        cObj[i] = -ru[1]*deltaPose[i];
    }
    // final velocity
    for( int i=0; i<NV; i++ )
    {
        cObj[6+i] = -ru[2]*deltaVelo[i];
    }
    // virtual stiffness
    dKConSN = 0;
    for( int i=0; i<NTS; i++ )
    {
        dKCon[i].setZero();
        for( int j=0; j<NPAIR; j++ )
        {
            dKCon[i][j] = -U[i][NU+j];
            cObj[6+NV+i*NPAIR+j] = -ru[3]*dKCon[i][j];
        }
//        cObj[6+i] = -ru[2]*dKCon[i].squaredNorm();
        dKConSN += dKCon[i].squaredNorm();
    }
    // constant objective term
    ObjAdd = 0.5*(ru[0]*deltaPose.block<2,1>(0,0).squaredNorm() +
                  ru[1]*deltaPose.block<4,1>(2,0).squaredNorm() +
                  ru[2]*deltaVelo.squaredNorm() +
                  ru[3]*dKConSN);
}

// solveCvx: solves the convex subproblem
void CitoSQOPT::solveCvx(double *xTraj, double r, const stateTraj X, const ctrlTraj U,
                         const stateDerTraj Fx, const ctrlDerTraj Fu, int *isJFree, int *isAFree,
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
    cvxProb.solve(Cold, qpHx, nc, n, neA, lencObj, nnH, iObj,
                  ObjAdd, valA, indA, locA, bl, bu, cObj,
                  eType, hs, x, pi, rc, nS, nInf, sInf, objective);
    // sort x to regular form
    this->sortX(x, indMove);
    // set only the trajectory-related variables (X and U)
    for( int i=0; i<NTRAJ; i++ )
    {
        xTraj[i] = x[i];
    }
    // cleaning
    delete []indA;  delete []locA; delete []valA;
    delete []x;     delete []bl;   delete []bu;
    delete []pi;    delete []rc;   delete []hs;   delete[]eType;
}

// setBounds: sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
void CitoSQOPT::setBounds(double r, const stateTraj X, const ctrlTraj U,
                          double *bl, double *bu, int *isJFree, int *isAFree,
                          double *qposLB, double *qposUB, double *tauLB, double *tauUB)
{
    // decision variables
    // * states
    for( int i=0; i<NTS+1; i++ )
    {
        for( int j=0; j<NV; j++ )
        {
            // ** change in free joint positions: unbounded
            // ** change in joint positions
            if( isJFree[j] == 0 )
            {
                bl[i*N+j] = qposLB[j] - X[i][j];
                bu[i*N+j] = qposUB[j] - X[i][j];
            }
            // ** change in joint velocities: unbounded (already set)
        }
        // * controls
        if( i < NTS )
        {
            // ** change in joint torques
            for( int j=0; j<NU; j++ )
            {
                if( isAFree == 0 )
                {
                    bl[dUOffset+i*M+j] = tauLB[j] - U[i][j];
                    bu[dUOffset+i*M+j] = tauUB[j] - U[i][j];
                }
            }
            // ** change in virtual stiffness
            for( int j=0; j< NPAIR; j++ )
            {
                bl[dUOffset+i*M+NU+j] = 0 - U[i][NU+j];
                bu[dUOffset+i*M+NU+j] = kCon0 - U[i][NU+j];
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
    for( int i=n; i<n+(NTS+1)*N; i++ )
    {
        bl[i] = 0;
        bu[i] = 0;
    }
    // 0 <= absolute value constraints <= inf
    for( int i=n+(NTS+1)*N; i<n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2; i++ )
    {
        bl[i] = 0;
    }
    // 0 <= trust region <= r
    bl[n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2] = 0;
    bu[n+(NTS+1)*N+((NTS+1)*N+NTS*M)*2] = r;
}

// setA: creates the sparse A matrix for linearized dynamics, auxiliary
// variables, and trust region constraints
void CitoSQOPT::setA(double *valA, int *indA, int *locA,
                     const stateDerTraj Fx, const ctrlDerTraj Fu)
{
    int colNo = 0, indNo = 0, indTS = 0;

    locA[0] = 0;
    // columns associated with dx[1,...,NTS]
    for( int i=0; i<NTS*N; i++ )
    {
        indA[indNo] = i;
        valA[indNo] = -1;
        indNo++;
        // time step index
        indTS = (int) floor(i/N);
        for( int j=0; j<N; j++ )
        {
            indA[indNo] = (indTS+1)*N+j;
            valA[indNo] = Fx[indTS](j,i%N);
            indNo++;
        }
        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;

        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with dx[NTS+1]
    for( int i=0; i<N; i++ )
    {
        indA[indNo] = NTS*N+i;
        valA[indNo] = -1;
        indNo++;

        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;

        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with du[1,...,NTS]
    for( int i=0; i<NTS*M; i++ )
    {
        // time step index
        indTS = (int) floor(i/M);
        for( int j=0; j<N; j++ )
        {
            indA[indNo] = (indTS+1)*N+j;
            valA[indNo] = Fu[indTS](j,i%M);
            indNo++;
        }

        // for auxiliary variables
        indA[indNo] = (NTS+1)*N + colNo;
        valA[indNo] = +1.0;
        indNo++; //auxNo1++;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + colNo;
        valA[indNo] = -1.0;
        indNo++; //auxNo2++;

        colNo++;
        locA[colNo] = indNo;
    }
    // columns associated with auxiliary variables
    int auxNo1 = 0, auxNo2 = 0;
    for( int i=0; i<(NTS+1)*N+NTS*M; i++ )
    {
        indA[indNo] = (NTS+1)*N + auxNo1;
        valA[indNo] = +1.0;
        indNo++; auxNo1++;
        indA[indNo] = (NTS+1)*N + (NTS+1)*N + NTS*M + auxNo2;
        valA[indNo] = +1.0;
        indNo++; auxNo2++;
        // for l1-norm
        indA[indNo] = (NTS+1)*N + 2*((NTS+1)*N + NTS*M);
        valA[indNo] = +1.0;
        indNo++;

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