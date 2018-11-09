// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //
#include "numdiff.h"
#include "scvx.h"
#include "savelog.h"
#include "snoptProblem.hpp"

// ******** quadratic cost H*x ************************************************/
void qpHx(int *nnH, double x[], double Hx[], int *nState,
    	       char   cu[], int   *lencu,
    	       int    iu[], int   *leniu,
    	       double ru[], int   *lenru)
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
    for( int i=6; i<6+NTS; i++ )
    {
      Hx[i] = ru[2]*x[i];
    }
}
// ********* create state/control vectors/matrices ****************************/
stateVecThread X, dX, XL;   ctrlVecThread  U, dU, Utemp;
stateMatThread Fx;          ctrlMatThread  Fu;
// ********* create objects ***************************************************/
SCvx sc;
NumDiff nd;
SaveLog sl;               FILE* logfile;
sqoptProblem cvxProb;
// ********* position & torque limits *****************************************/
double *qpos_lb = new double[NV];     double *qpos_ub = new double[NV];
double *tau_lb  = new double[ndof];   double *tau_ub  = new double[ndof];
int    *isObj   = new int[NV];
// ********* extra variables **************************************************/
Eigen::VectorXd po_d(6), dpo_d(6);
Eigen::VectorXd dkcon_d(NTS), dvel_d((NTS+1)*NV);
// ********* sqopt initialization *********************************************/
int nnH     = 6 + NTS;
int lencObj = 6 + NTS;
int lenru = 4;
double *cObj = new double[lencObj]; double *ru   = new double[lenru];
int neA = NTS*N*N + (NTS+1)*N + NTS*N*M + ((NTS+1)*N+NTS*M)*5;
int n   = ((NTS+1)*N + NTS*M)*2;              // x2 is for auxiliary variables
int nc  = (NTS+1)*N + ((NTS+1)*N+NTS*M)*2 + 1;
int iObj = -1;
int nS, nInf;
double ObjAdd  = 0, sInf = 0, objective;
double infBnd = 1.0e20;
int Cold = 0, Basis = 1, Warm = 2;
// ********* mujoco initialization ********************************************/
mjModel*  m = NULL;
mjData*   d = NULL;
//==============================================================================
int main(int argc, char const *argv[]) {
    // activate mujoco
    mj_activate("/home/aykut/Development/cito/src/mjc200/bin/mjkey.txt");
    // load xml model
    m = mj_loadXML("/home/aykut/Development/cito/src/mjc200/model/flymanoid.xml", 0, 0, 0);
    if( !m )
        mju_error("Cannot load the model");
    // initial pose
    qpos0 << -M_PI/12, -M_PI/beta, 0, 3*M_PI/beta, -2*M_PI/12, -2*M_PI/beta, 0;
    // get position & torque limits
    sc.getBounds(qpos_lb, qpos_ub, tau_lb, tau_ub, isObj, m);
    // allocate state and control matrices/vectors
    X.resize(NTS+1);  U.resize(NTS);  Utemp.resize(NTS);
    Fx.resize(NTS);   Fu.resize(NTS);
    dX.resize(NTS+1); dU.resize(NTS); XL.resize(NTS+1);

    stateMatThread Fxtest; Fxtest.resize(NTS);
    ctrlMatThread  Futest; Futest.resize(NTS);
    // initialize & set options for sqopt
    cvxProb.initialize("", 1);
    // cvxProb.setSpecsFile("/home/aykut/Development/cito_manipulation/src/mujoco/src/sqoptOptions.spc");
    // ********* objective ****************************************************/
    // desired position

    // FIXME FIXME FIXME
    // desired z needs to be removed!
    // FIXME FIXME FIXME
    po_d << 1.25, 0, 0.0, 0.0, 0.0, 0.0; dpo_d.setZero();

    for( int i=0; i<2; i++ )
    {
      dpo_d[i] = po_d[i] - m->qpos0[i];
    }
    // weights
    ru[0] = 1e3;
    ru[1] = 1e0;
    ru[2] = 2e-2;
    cvxProb.setUserR(ru, lenru);
    // linear and constant terms
    // Hessian indices
    int nMove = nnH;
    int *indMove = new int[nMove];
    // k variables
    for( int i=0; i < NTS; i ++ )
    {
      indMove[i] = (NTS+1)*N + NTS*M - 1 - (i*M);
    }
    for(int i=0; i<6; i++ )
    {
        indMove[nnH-i-1] = NTS*N+i;
    }
    // for(int i=0;i<nnH; i++)
    // std::cout << "indMove" << i << ": " << indMove[i] << '\n';
    // ********* initial control guess ****************************************/
    ctrlVec_t U0; U0.setZero();
    for( int i=0; i<NTS; i++ )
    {
      U[i] = U0;
      for( int j=0; j<ncc; j++ )
      {   U[i](ndof) = kcon0; }
    }
    // ********* scvx parameters **********************************************/
    int maxIter = 25;
    double Jtemp[maxIter+1], J[maxIter+1], L[maxIter+1];
    double r[maxIter+1], rho[maxIter+1], dL[maxIter+1], dJ[maxIter+1];
    r[0] = 1e2;   // initial trust region radius
    double dLTol = 1e-4;
    double rho0 = 0, rho1 = 0.25, rho2 = 0.90, rMin = 0, rMax = 1e20;
    double alpha = 2, beta = 3.2;
    int accept[maxIter+1], dLTolMet = 0, stop = 0;
    // scvx ====================================================================
    int iter = 0;
    while( stop == 0 )
    {
        // clean state and linearization matrices
        for( int i=0; i<NTS; i++ )
        {
          X[i].setZero();   XL[i].setZero();  dX[i].setZero();
          Fx[i].setZero();  Fu[i].setZero();  dU[i].setZero(); Utemp[i].setZero();
        }
        X[NTS].setZero();   XL[NTS].setZero();  dX[NTS].setZero();
        // simulation ==========================================================
        // make mjData
        d = mj_makeData(m);
        // set initial pose
        mju_copy(d->qpos+nobj*7, qpos0.data(), ndof);
        mj_forward(m, d);
        sc.setControl(m, d, U[0]);
        // rollout and linearize the dynamics
        for( int i=0; i<NTS; i++ )
        {
            // get the current state values
            X[i] = nd.getState(m, d);
            //linearization
            nd.cvxDyn(Fx[i].data(), Fu[i].data(), m, d, U[i]);
            // take tc/dt steps
            for( int j=0; j<ndpc; j++ )
            {
              // initialize the step
              mj_step1(m, d);
              // set ctrl and xfrc
              sc.setControl(m, d, U[i]);
              // complete the step
              mj_step2(m, d);
          }
        }
        X[NTS] = nd.getState(m, d);
        // delete data
        mj_deleteData(d);
        // get the nonlinear cost if first iteration
        if( iter == 0 ) { J[iter] = sc.getCost(X, U, po_d, ru); }
        // optimization ========================================================
        // fresh start
        int    *indA = new int[neA];
        int    *locA = new int[n+1];
        double *valA = new double[neA];
        double *x  = new double[n+nc];
        double *bl = new double[n+nc];
        double *bu = new double[n+nc];
        double *pi = new double[nc];
        double *rc = new double[n+nc];
        int    *hs = new    int[n+nc];
        int *eType = new    int[n+nc];
        // *********** initial guess ******************************************/
        for( int i=0; i<n+nc; i++ )
        {
          x[i]  = 0.0;
          bl[i] = -infBnd;  bu[i] = +infBnd;
          hs[i] = 0;  eType[i] = 0;  rc[i] = 0.0;
          if( i>=n ) { pi[i-n] = 0.0; }
        }
        // *********** set linear constraints and bounds **********************/
        sc.setA(valA, indA, locA, Fx, Fu);
        sc.setBounds(bl, bu, qpos_lb, qpos_ub, tau_lb, tau_ub, isObj, n, nc, X, U, r[iter]);
        int less_counter = 0, iMove = 0;
        for( int i=0; i<nnH; i++ )
        {
          if( i>0 && indMove[i]<indMove[i-1] ) { less_counter += 1; }
          iMove = indMove[i]+less_counter;

          sc.moveColA(valA, indA, locA, iMove, neA, n);
          sc.moveRowBounds(bl, bu, iMove);
        }
        // *********** set weights ********************************************/
        // ru[0] = ru[0];
        // cvxProb.setUserR(ru, lenru);
        // *********** set linear and constant objective terms ****************/
        for( int i=0; i<6; i++ )
        {
          dpo_d[i] = po_d[i] - X[NTS](i);
        }
        for( int i=0; i<2; i++ )
        {
          cObj[i] = -ru[0]*dpo_d[i];
        }
        for( int i=2; i<6; i++ )
        {
          cObj[i] = -ru[1]*dpo_d[i];
        }
        // kcon terms
        for( int i=0; i<NTS; i++ )
        {
          dkcon_d[i] = -U[i](ndof);
          cObj[6+i] = -ru[2]*dkcon_d[i];
        }
        ObjAdd = 0.5*(ru[0]*dpo_d.block<2,1>(0,0).squaredNorm() +
                      ru[1]*dpo_d.block<4,1>(2,0).squaredNorm() +
                      ru[2]*dkcon_d.squaredNorm());
        // *********** sqopt solve ********************************************/
        cvxProb.solve(Cold, qpHx, nc, n, neA, lencObj, nnH, iObj,
                      ObjAdd, valA, indA, locA, bl, bu, cObj,
              	      eType, hs, x, pi, rc, nS, nInf, sInf, objective);
        // *********** sort solution*******************************************/
        sc.sortX(x, indMove, nMove, n);
        // *********** linear cost ********************************************/
        for( int i=0; i<NTS+1; i++ )
        {
          for( int j=0; j<N; j++ )
          {
            dX[i](j) = x[i*N+j];
          }
          if( i < NTS )
          {
            for( int j=0; j<M; j++ )
            {
              dU[i](j) = x[(NTS+1)*N+i*M+j];
            }
          }
          XL[i] = X[i] + dX[i];
        }
        // *********** nonlinear cost *****************************************/
        // make mjData
        d = mj_makeData(m);
        // set initial pose
        mju_copy(d->qpos+nobj*7, qpos0.data(), ndof);
        mj_forward(m, d);
        sc.setControl(m, d, U[0]);
        // rollout the dynamics with the optimal change in controls
        for( int i=0; i<NTS; i++ )
        {
            // get the current state values
            X[i] = nd.getState(m, d);
            // temporary new control
            Utemp[i] = U[i]+dU[i];
            // take tc/dt dynamic steps
            for( int j=0; j<ndpc; j++ )
            {
              // initialize the step
              mj_step1(m, d);
              // set ctrl and xfrc
              sc.setControl(m, d, Utemp[i]);
              // complete the step
              mj_step2(m, d);
          }
        }
        X[NTS] = nd.getState(m, d);
        std::cout << "\n\nIteration " << iter << ":" << '\n';
        std::cout << "X:  " << X[NTS].transpose() << "\n";
        std::cout << "XL: " << XL[NTS].transpose() << "\n";
        std::cout << "dX: " << dX[NTS].transpose() << "\n\n";
        // get the linear and nonlinear costs
        L[iter]     = sc.getCost(XL, Utemp, po_d, ru);
        Jtemp[iter] = sc.getCost(X,  Utemp, po_d, ru);
        // *********** similarity measure *************************************/
        dJ[iter] = J[iter] - Jtemp[iter];
        dL[iter] = J[iter] - L[iter];
        rho[iter] = dJ[iter]/dL[iter];
        if( dL[iter] > 0 && dL[iter] < dLTol )
        {
          dLTolMet = 1;
        }
        // *********** accept/reject solution *********************************/
        if( rho[iter]<=rho0 || (dL[iter]<0 && dJ[iter]<0) || std::isnan(rho[iter]) )
        {
          accept[iter] = 0;
          r[iter+1] = r[iter]/alpha;
          J[iter+1] = J[iter];
        }
        else { accept[iter] = 1; }
        if( accept[iter] )
        {
          J[iter+1] = Jtemp[iter];
          for( int i=0; i<NTS; i++ )
          {
            U[i] = Utemp[i];
          }
          if( rho[iter] < rho1 )
          { r[iter+1] = r[iter]/alpha;  }
          else if( rho[iter]>=rho1 && rho[iter]<rho2 )
          { r[iter+1] = r[iter];        }
          else if( rho[iter]>=rho2 )
          { r[iter+1] = r[iter]*beta;   }
        }
        // bound r
        r[iter+1] = std::max(r[iter+1], rMin);  // lower bound
        r[iter+1] = std::min(r[iter+1], rMax);  // upper bound
        // *********** temporary clean up *************************************/
        // sqopt
        delete []indA;  delete []locA; delete []valA;
        delete []x;     delete []bl;   delete []bu;
        delete []pi;    delete []rc;   delete []hs;   delete[]eType;
        // delete data
        mj_deleteData(d);
        // *********** next iteration *****************************************/
        iter += 1;
        // *********** stop condition *****************************************/
        if( iter == maxIter )
        {
          stop = 1;
          std::cout << "\n\n\tWARNING: Maximum number of iterations reached\n\n";
        }
        if( dLTolMet!=0 )
        {
          stop = 1;
          std::cout << "\n\n\tWARNING: dL<dLTol\n\n";
        }
    }
    // *********** print iteration details ************************************/
    std::cout << "\n\nSCVX Summary\nJ0=" << J[0] << "\n\n";
    for( int i=0; i<iter; i++ )
    {
      if( i%10 == 0 )
      {
        printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
        "Iteration","L","J","dL","dJ","rho","r","accept");
      }
      printf("%-12d%-12.6g%-12.6g%-12.3g%-12.3g%-12.3g%-12.3g%-12d\n",
      i,L[i],Jtemp[i],dL[i],dJ[i],rho[i],r[i],accept[i]);
    }
    // *********** record resulting motion ************************************/
    // make mjData
    d = mj_makeData(m);
    // set initial pose
    mju_copy(d->qpos+nobj*7, qpos0.data(), ndof);
    mj_forward(m, d);
    sc.setControl(m, d, U[0]);
    // open log file and write header
    logfile = fopen("recorded_traj_flymanoid_scvx", "wb");
    if (logfile == NULL) { mju_error("ERROR:  Unable to open file"); }
    sl.writeHeader(m, d, logfile);
    // rollout the dynamics
    for( int i=0; i<NTS; i++ )
    {
        // get the current state values
        X[i] = nd.getState(m, d);
        std::cout << "U" << i << " = " << U[i].transpose() << '\n';
        // take tc/dt dynamic steps
        for( int j=0; j<ndpc; j++ )
        {
          // initialize the step
          mj_step1(m, d);
          // set ctrl and xfrc
          sc.setControl(m, d, U[i]);
          // complete the step
          mj_step2(m, d);
          sl.writeData(m, d, logfile);
      }
    }
    X[NTS] = nd.getState(m, d);
    std::cout << "\nFinal state:" << '\n';
    std::cout << "X:  " << X[NTS].transpose() << "\n";
    std::cout << "XL: " << XL[NTS].transpose() << "\n";
    std::cout << "dX: " << dX[NTS].transpose() << "\n\n";
    // delete data and close log file
    mj_deleteData(d);
    fclose(logfile);
    // *********** shut down **************************************************/
    // custom variables
    delete []qpos_lb; delete []qpos_ub; delete [] isObj;
    delete []tau_lb;  delete []tau_ub;  delete [] ru;
    // mujoco
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
