// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SCVX class consists of functions that setup the successive
// convexification algorithm.

#include "cito_scvx.h"

// ***** CONSTRUCTOR ***********************************************************
CitoSCvx::CitoSCvx(const mjModel* model) : m(model), cc(model), nd(model)
{
    r[0] = r0;
    X.resize(NTS+1);    dX.resize(NTS+1);   XL.resize(NTS+1);
    U.resize(NTS);      dU.resize(NTS);     Utemp.resize(NTS);
    Fx.resize(NTS);     Fu.resize(NTS);
}

// getCost: returns the nonlinear cost given state matrix X
double CitoSCvx::getCost(const stateVecThread X, const ctrlVecThread U, const Eigen::VectorXd po_d, const double *w)
{
    // terminal cost
    for( int i=0; i<6; i++ )
    {
        po_f[i] = X[NTS](i);
    }
    Jt = 0.5*(w[0]*(po_d.block<2,1>(0,0)-po_f.block<2,1>(0,0)).squaredNorm()+
              w[1]*(po_d.block<4,1>(2,0)-po_f.block<4,1>(2,0)).squaredNorm());
    // integrated cost
    k.setZero();
    for( int i=0; i<NTS; i++ )
    {
        for( int j=0; j<params::npair; j++ )
        {
            k[i,j] = U[i](j);
        }
    }
    Ji = 0.5*w[2]*k.squaredNorm();
    // total cost
    J = Jt + Ji;

    return J;
}

// runSimulation: runs forward simulation given a control trajectory
void CitoControl::runSimulation(const ctrlMatThread U)
{
    // clean state and linearization matrices
    for( int i=0; i<NTS; i++ )
    {
        X[i].setZero();   XL[i].setZero();  dX[i].setZero();
        Fx[i].setZero();  Fu[i].setZero();  dU[i].setZero(); Utemp[i].setZero();
    }
    X[NTS].setZero();   XL[NTS].setZero();  dX[NTS].setZero();
    // make mjData
    d = mj_makeData(m);
    // initialize d
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    cc.setControl(m, d, U[0]);
    // rollout and linearize the dynamics
    for( int i=0; i<NTS; i++ )
    {
        // get the current state values
        X[i] = this->getState(m, d);
        // linearization
        nd.linDyn(Fx[i].data(), Fu[i].data(), m, d, U[i]);
        // take tc/dt steps
        cc.takeStep(d, U[i]);
    }
    X[NTS] = nd.getState(m, d);
    // delete data
    mj_deleteData(d);
}

void CitoSCvx::solve()
{
    int iter = 0;
    while( stop == 0 )
    {
        // simulation ==========================================================

        // get the nonlinear cost if first iteration
        if( iter == 0 ) { J[iter] = this->getCost(X, U, po_d, ru); }
        // optimization ========================================================
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
        // *********** initial guess ******************************************/
        for( int i=0; i<n+nc; i++ )
        {
            x[i]  = 0.0;
            bl[i] = -infBnd;  bu[i] = +infBnd;
            hs[i] = 0;  eType[i] = 0;  rc[i] = 0.0;
            if( i>=n ) { pi[i-n] = 0.0; }
        }
        // *********** set linear constraints and bounds **********************/
        sq.setA(valA, indA, locA, Fx, Fu);
        sq.setBounds(bl, bu, cc.qpos_lb, cc.qpos_ub, cc.tau_lb, cc.tau_ub, cc.isJFree, cc.isAFree, n, nc, X, U, r[iter]);

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
        X[NTS] = cc.getState(m, d);
        std::cout << "\n\nIteration " << iter << ":" << '\n';
        std::cout << "X:  " << X[NTS].transpose() << "\n";
        std::cout << "XL: " << XL[NTS].transpose() << "\n";
        std::cout << "dX: " << dX[NTS].transpose() << "\n\n";
        // get the linear and nonlinear costs
        L[iter]     = this->getCost(XL, Utemp, po_d, ru);
        Jtemp[iter] = this->getCost(X,  Utemp, po_d, ru);
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
}