// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SQOPT class consists of functions that setup the CITO problem for the
// convex programming solver SQOPT.

// ***** CLASS TYPE ************************************************************
// Solver specific

#include "cito_sqopt.h"

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoSQOPT::CitoSQOPT()
{
    dKcon.resize(NTS);
    // indices to move
    // *** virtual stiffness variables
    for( int i=0; i < NTS; i ++ )
    {
        for( int j=0; j<NPAIR; j++ )
        {
            indMove[i*NPAIR+j] = (NTS+1)*N + NTS*M - 1 - (i*M + j);
        }
    }
    // *** pose variables for the control joint
    for( int i=0; i<6; i++ )
    {
        indMove[nnH-1-i] = NTS*N + task::controlJointPos0 + i;
    }
    // initialize & set options for SQOPT
    cvxProb.initialize("", 1);
    // set the weights
    ru[0] = task::w1; ru[1] = task::w2; ru[2] = task::w3;
    cvxProb.setUserR(ru, lenru);
}
// ***** FUNCTIONS *************************************************************
// qpHx: sets Hx to the H*x part of the quadratic cost to be multiplied by x'
void qpHx(int *nnH, double x[], double Hx[], int *nState,
          char cu[], int *lencu, int iu[], int *leniu,
          double ru[], int *lenru)
{
    // final position of object
    for( int i=0; i<2; i++ )
    {
        Hx[i] = ru[0]*x[i];
    }
    for( int i=2; i<6; i++ )
    {
        Hx[i] = ru[1]*x[i];
    }
    // kcon terms
    for( int i=6; i<6+NTS*NPAIR; i++ )
    {
        Hx[i] = ru[2]*x[i];
    }
}

// solveCvx: solves the convex subproblem
void CitoSQOPT::solveCvx(double *xTraj, double r, const stateVecThread X, const ctrlVecThread U,
                         const stateDerThread Fx, const ctrlDerThread Fu, int *isJFree, int *isAFree,
                         double *qpos_lb, double *qpos_ub, double *tau_lb, double *tau_ub)
{
    // fresh start
    int    *indA  = new int[neA];
    int    *locA  = new int[n+1];
    double *valA  = new double[neA];
    double *x     = new double[n+nc];
    double *bl    = new double[n+nc];
    double *bu    = new double[n+nc];
    double *pi    = new double[nc];
    double *rc    = new double[n+nc];
    int    *hs    = new    int[n+nc];
    int    *eType = new    int[n+nc];
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
    this->setBounds(r, X, U, bl, bu, isJFree, isAFree, qpos_lb, qpos_ub, tau_lb, tau_ub);
    // sort A and bounds w.r.t. the order in qpHX
    this->sortToMatch(valA, indA, locA, indMove, bl, bu);
    // set linear and constant cost terms
    this->setCost(X, U, ru, cObj, ObjAdd);
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

// setCost: sets linear and constant cost terms
void CitoSQOPT::setCost(const stateVecThread X, const ctrlVecThread U,
                        double *ru, double *cObj, double& ObjAdd)
{
    // set weights
     cvxProb.setUserR(ru, lenru);
    // *********** set linear and constant objective terms ****************/
    // desired change in the pose
    dPose.setZero();
    for( int i=0; i<6; i++ )
    {
        dPose[i] = task::desiredPose[i] - X[NTS][task::controlJointPos0+i];
    }
    // position
    for( int i=0; i<2; i++ )
    {
        cObj[i] = -ru[0]*dPose[i];
    }
    // orientation
    for( int i=2; i<6; i++ )
    {
        cObj[i] = -ru[1]*dPose[i];
    }
    // virtual stiffness
    dKconSN = 0;
    for( int i=0; i<NTS; i++ )
    {
        dKcon[i].setZero();
        for( int j=0; j<NPAIR; j++ )
        {
            dKcon[i][j] = -U[i][NU+j];
        }
        cObj[6+i] = -ru[2]*dKcon[i].squaredNorm();
        dKconSN += dKcon[i].squaredNorm();
    }
    // constant objective term
    ObjAdd = 0.5*(ru[0]*dPose.block<2,1>(0,0).squaredNorm() +
                  ru[1]*dPose.block<4,1>(2,0).squaredNorm() +
                  ru[2]*dKconSN);
}

// setBounds: sets bounds of dX, dU, and constraints (dynamics, trust region, etc.)
void CitoSQOPT::setBounds(double r, const stateVecThread X, const ctrlVecThread U,
                          double *bl, double *bu, int *isJFree, int *isAFree,
                          double *qpos_lb, double *qpos_ub, double *tau_lb, double *tau_ub)
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
                bl[i*N+j] = qpos_lb[j] - X[i](j);
                bu[i*N+j] = qpos_ub[j] - X[i](j);
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
                    bl[dU_offset+i*M+j] = tau_lb[j] - U[i](j);
                    bu[dU_offset+i*M+j] = tau_ub[j] - U[i](j);
                }
            }
            // ** change in virtual stiffness
            for( int j=0; j< NPAIR; j++ )
            {
                bl[dU_offset+i*M+NU+j] = 0 - U[i](NU+j);
                bu[dU_offset+i*M+NU+j] = params::kcon0[j] - U[i](NU+j);
            }
        }
    }
    // auxiliary variables > 0
    for( int i=aux_offset; i<n; i++ )
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
                     const stateDerThread Fx, const ctrlDerThread Fu)
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
        indTS = floor(i/N);
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
        indTS = floor(i/M);
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

// moveColA: moves iMove to left
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

// moveRowBounds: moves iMove to top
void CitoSQOPT::moveRowBounds(double *bl, double *bu, int iMove)
{
    double bl_temp[iMove], bu_temp[iMove];
    bl_temp[0] = bl[iMove];
    bu_temp[0] = bu[iMove];
    for( int i=0; i<iMove+1; i++ )
    {
        if( i<iMove )
        {
            bl_temp[i+1] = bl[i];
            bu_temp[i+1] = bu[i];
        }
        bl[i] = bl_temp[i];
        bu[i] = bu_temp[i];
    }
}

// sortX: sorts decision variables back to original
void CitoSQOPT::sortX(double *x, int *moveIndices)
{
    int k = 0;
    for( int i=0; i<n; i++ )
    {
        xtemp[i] = x[nnH+k];
        for( int j=0; j<nnH; j++ )
        {
            if( i == moveIndices[j] )
            {
                xtemp[moveIndices[j]] = x[nnH-j-1];
                k = k - 1;
            }
        }
        k = k + 1;
    }
    for( int i=0; i<n; i++ )
    {
        x[i] = xtemp[i];
    }
}

// Test
// // print A
// int neCol[n];
// for( int i=0; i<n; i++ )
// {
//   neCol[i] = locA[i+1]-locA[i];
// }
// int indNo = 0, colNo = 0;
// for( int i=0; i<n; i++ )
// {
//   std::cout << "Column " << i << ":\n";
//   for( int j=0; j<neCol[i]; j++ )
//   {
//     std::cout << "\tindNo: " << indNo << ", indA: " << indA[indNo] << ", valA: " << valA[indNo] << '\n';
//     indNo+=1;
//   }
//   colNo++;
//   std::cout << "\tlocA: " << locA[colNo] << '\n';
// }
// std::cout << "n: " << n << ", neA: " << neA << '\n';

// // test shifts
// int ntest = 10;
// double *test   = new double[ntest];
// double *testDummy  = new double[ntest];
// double assignTest[10] = {1, 2, 3, 4, 5, 6 ,7, 8, 9, 10};
// for( int i=0; i<ntest; i++ )
// {
//   test[i] = assignTest[i];
//   testDummy[i] = assignTest[i];
// }
// for( int i=0; i<ntest; i++ )
//   std::cout << test[i] << ' ';
// std::cout << "\n\n\n";
// const int nshift = 5;
// int *shift = new int[nshift];
// int assignShift[nshift] = {3, 9, 4, 7, 8};
// for( int i=0; i<nshift; i++ )
//   shift[i] = assignShift[i]-1;
// int less_counter = 0;
// for( int i=0; i<nshift; i++ )
// {
//   if( i>0 && shift[i]<shift[i-1] )
//   {
//     less_counter++;
//   }
//   int rMove = shift[i]+less_counter;
//   sc.moveRowBounds(test, testDummy, rMove);
//   for( int i=0; i<ntest; i++ )
//     std::cout << test[i] << ' ';
//   std::cout << "\n\n\n";
// }
//
// sc.sortX(test, shift, nshift, ntest);
// for( int i=0; i<ntest; i++ )
//   std::cout << test[i] << ' ';
// std::cout << "\n\n\n";