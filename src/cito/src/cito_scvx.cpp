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
    Fx.resize(NTS);     Fu.resize(NTS);     Kcon.resize(NTS);
}

// getCost: returns the nonlinear cost given state matrix X
double CitoSCvx::getCost(const stateVecThread X, const ctrlVecThread U, const Eigen::VectorXd linFlag, const double *w)
{
    // terminal cost
    for( int i=0; i<6; i++ )
    {
        finalPose[i] = X[NTS](i);
    }
    Jt = 0.5*(w[0]*(desiredPose.block<2,1>(0,0)-finalPose.block<2,1>(0,0)).squaredNorm()+
              w[1]*(desiredPose.block<4,1>(2,0)-finalPose.block<4,1>(2,0)).squaredNorm());
    // integrated cost
    KconSN = 0;
    for( int i=0; i<NTS; i++ )
    {
        Kcon[i].setZero();
        for( int j=0; j<NPAIR; j++ )
        {
            Kcon[i][j] = U[i][NU+j];
        }
        KconSN += Kcon[i].squaredNorm();
    }
    Ji = 0.5*w[2]*KconSN;
    // total cost
    J = Jt + Ji;

    return J;
}

// runSimulation: rollouts and linearizes the dynamics given a control trajectory
void CitoSCvx::runSimulation(const ctrlMatThread U, bool linearize)
{
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
        X[i].setZero();
        X[i] = cc.getState(d);
        // linearization
        if( linearize )
        {
            Fx[i].setZero(); Fu[i].setZero();
            nd.linDyn(d, U[i], Fx[i].data(), Fu[i].data());
        }
        // take tc/dt steps
        cc.takeStep(d, U[i]);
    }
    X[NTS].setZero();
    X[NTS] = cc.getState(d);
    // delete data
    mj_deleteData(d);
}

void CitoSCvx::solve(const ctrlMatThread U)
{
    int iter = 0;
    while( stop == 0 )
    {
        // simulation ==========================================================
        this->runSimulation(U, true);
        // get the nonlinear cost if the first iteration
        if( iter == 0 ) { J[iter] = this->getCost(X, U, po_d, ru); }
        // convex optimization =================================================
        // solve for the optimal change in trajectory
        double *dTraj = new double[NTRAJ];
        sq.solveCVX(dTraj, r[iter], X, U, Fx, Fu, qpos_lb, qpos_ub, tau_lb, tau_ub, isJFree, isAFree);
        // apply the change
        for( int i=0; i<NTS+1; i++ )
        {
            dX[i].setZero(); XL[i].setZero();
            for( int j=0; j<N; j++ )
            {
                dX[i][j] = dTraj[i*N+j];
            }
            if( i < NTS )
            {
                dU[i].setZero(); Utemp[i].setZero();
                for( int j=0; j<M; j++ )
                {
                    dU[i][j] = dTraj[(NTS+1)*N+i*M+j];
                }
                Utemp[i] = U[i] + dU[i];
            }
            XL[i] = X[i] + dX[i];
        }
        // evaluate the dynamics for the change and get the cost values ========
        this->runSimulation(Utemp, false);
        // get the linear and nonlinear costs
        L[iter]     = this->getCost(XL, Utemp, desiredPose, ru);
        Jtemp[iter] = this->getCost(X,  Utemp, desiredPose, ru);
        // similarity measure ==================================================
        dJ[iter] = J[iter] - Jtemp[iter];
        dL[iter] = J[iter] - L[iter];
        rho[iter] = dJ[iter]/dL[iter];
        if( dL[iter] > 0 && dL[iter] < dLTol )
        {
            dLTolMet = 1;
        }
        // accept or reject the solution =======================================
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
        // next iteration ======================================================
        iter++;
        // stopping criteria check =============================================
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
        // screen output for the iteration =====================================
        std::cout << "\n\nIteration " << iter << ":" << '\n';
        std::cout << "X:  " << X[NTS].transpose() << "\n";
        std::cout << "XL: " << XL[NTS].transpose() << "\n";
        std::cout << "dX: " << dX[NTS].transpose() << "\n\n";
    }
    // *********** screen output for the whole process ************************/
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