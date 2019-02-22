// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SCVX class consists of functions that setup and execute the successive
// convexification algorithm.

#include "cito_scvx.h"

// ***** CONSTRUCTOR ***********************************************************
CitoSCvx::CitoSCvx(const mjModel* model) : m(model), cc(model), nd(model)
{
    // read task parameters
    YAML::Node parameters = YAML::LoadFile(paths::taskConfig);
    std::vector<double> desiredPoseInput = { parameters["desiredPoseInput"].as<std::vector<double>>() };
    std::vector<double> desiredVeloInput = { parameters["desiredVeloInput"].as<std::vector<double>>() };
    desiredPose = Eigen::Map<Eigen::Matrix<double, 6, 1>>(desiredPoseInput.data(), desiredPoseInput.size());
    desiredVelo = Eigen::Map<Eigen::Matrix<double, 6, 1>>(desiredVeloInput.data(), desiredVeloInput.size());
    controlJointDOF0 = parameters["controlJointDOF0"].as<int>();
    // initial trust region radius
    r[0] = r0;
    // set bounds
    cc.getBounds();
    // trajectories
    XSucc.resize(NTS+1);    dX.resize(NTS+1);   XTilde.resize(NTS+1);
    USucc.resize(NTS);      UTemp.resize(NTS);  dU.resize(NTS);
    Fx.resize(NTS);         Fu.resize(NTS);
    KCon.resize(NTS);
    // read task parameters
    weight[0] = parameters["w1"].as<double>();
    weight[1] = parameters["w2"].as<double>();
    weight[2] = parameters["w3"].as<double>();
    weight[3] = parameters["w4"].as<double>();
}

// ***** FUNCTIONS *************************************************************
// getCost: returns the nonlinear cost given state matrix X
double CitoSCvx::getCost(stateVec_t XFinal, const ctrlVecThread U)
{
    // terminal cost
    for( int i=0; i<6; i++ )
    {
        finalPose[i] = XFinal[controlJointDOF0 + i];
        finalVelo[i] = XFinal[controlJointDOF0 + NV + i];
    }
    Jf = 0.5*(weight[0]*(desiredPose.block<2,1>(0,0)-finalPose.block<2,1>(0,0)).squaredNorm()+
              weight[1]*(desiredPose.block<4,1>(2,0)-finalPose.block<4,1>(2,0)).squaredNorm()+
              weight[3]*(desiredVelo - finalVelo).squaredNorm());
    // integrated cost
    KConSN = 0;
    for( int i=0; i<NTS; i++ )
    {
        KCon[i].setZero();
        for( int j=0; j<NPAIR; j++ )
        {
            KCon[i][j] = U[i][NU+j];
        }
        KConSN += KCon[i].squaredNorm();
    }
    Ji = 0.5*weight[2]*KConSN;
    // total cost
    Jt = Jf + Ji;
    return Jt;
}

// runSimulation: rollouts and linearizes the dynamics given a control trajectory
trajectory CitoSCvx::runSimulation(const ctrlVecThread U, bool linearize, bool save)
{
    // make mjData
    mjData* d = NULL;
    d = mj_makeData(m);
    // initialize d
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    cc.setControl(d, U[0]);
    // rollout (and linearize) the dynamics
    for( int i=0; i<NTS; i++ )
    {
        // get the current state values
        XSucc[i].setZero();
        XSucc[i] = cc.getState(d);
        // linearization
        if( linearize )
        {
            Fx[i].setZero(); Fu[i].setZero();
            nd.linDyn(d, U[i], Fx[i].data(), Fu[i].data());
        }
        // take tc/dt steps
        cc.takeStep(d, U[i], save);
    }
    XSucc[NTS].setZero();
    XSucc[NTS] = cc.getState(d);
    // delete data
    mj_deleteData(d);
    // build trajectory
    traj.X = XSucc; traj.U = U;
    if( linearize )
    {
        traj.Fx = Fx; traj.Fu = Fu;
    }
    return traj;
}

// solveSCvx: executes the successive convexification algorithm
ctrlVecThread CitoSCvx::solveSCvx(const ctrlVecThread U0)
{
    // initialize USucc for the first succession
    for( int i=0; i<NTS; i++ ) { USucc[i] = U0[i]; }
    // start the SCvx algorithm
    int iter = 0;
    while( stop == 0 )
    {
        // simulation and linearization ========================================
        trajS = {};
        trajS = this->runSimulation(USucc, true, false);
        // get the nonlinear cost if the first iteration
        if( iter == 0 ) { J[iter] = this->getCost(trajS.X[NTS], USucc); }
        // convex optimization =================================================
        double *dTraj = new double[NTRAJ];
        sq.solveCvx(dTraj, r[iter], trajS.X, USucc, trajS.Fx, trajS.Fu, cc.isJFree, cc.isAFree,
                    cc.qposLB, cc.qposUB, cc.tauLB, cc.tauUB);
        // apply the change
        for( int i=0; i<NTS+1; i++ )
        {
            // states
            dX[i].setZero(); XTilde[i].setZero();
            for( int j=0; j<N; j++ )
            {
                dX[i][j] = dTraj[i*N+j];
            }
            XTilde[i] = trajS.X[i] + dX[i];
            // controls
            if( i < NTS )
            {
                dU[i].setZero(); UTemp[i].setZero();
                for( int j=0; j<M; j++ )
                {
                    dU[i][j] = dTraj[(NTS+1)*N+i*M+j];
                }
                UTemp[i] = USucc[i] + dU[i];
            }
        }
        // evaluate the dynamics for the change and get the cost values ========
        trajTemp = {};
        trajTemp = this->runSimulation(UTemp, false, false);
        // get the linear and nonlinear costs
        JTilde[iter]     = this->getCost(XTilde[NTS], UTemp);
        JTemp[iter] = this->getCost(trajTemp.X[NTS], UTemp);
        // similarity measure ==================================================
        dJ[iter] = J[iter] - JTemp[iter];
        dL[iter] = J[iter] - JTilde[iter];
        rho[iter] = dJ[iter]/dL[iter];
        if( dL[iter] > 0 && dL[iter] < dLTol )
        {
            dLTolMet = 1;
        }
        // accept or reject the solution =======================================
        // reject
        if( rho[iter]<=rho0 || (dL[iter]<0 && dJ[iter]<0) || std::isnan(rho[iter]) )
        {
            accept[iter] = 0;
            r[iter+1] = r[iter]/alpha;
            J[iter+1] = J[iter];
        }
        else { accept[iter] = 1; }
        // accept
        if( accept[iter] )
        {
            J[iter+1] = JTemp[iter];
            for( int i=0; i<NTS; i++ )
            {
                USucc[i] = UTemp[i];
            }
            if( rho[iter] < rho1 )
            { r[iter+1] = r[iter]/alpha;  }
            else if( rho[iter]>=rho1 && rho[iter]<rho2 )
            { r[iter+1] = r[iter];        }
            else if( rho[iter]>=rho2 )
            { r[iter+1] = r[iter]*beta;   }
        }
        // bound the trust region radius r
        r[iter+1] = std::max(r[iter+1], rMin);
        r[iter+1] = std::min(r[iter+1], rMax);
        // stopping criteria check =============================================
        if( iter+1 == maxIter )
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
        std::cout << "\n\nIteration " << iter+1 << ":" << '\n';
        std::cout << "X:      " << trajTemp.X[NTS].transpose() << "\n";
        std::cout << "XTilde: " << XTilde[NTS].transpose() << "\n";
        std::cout << "JTilde = " << JTilde[iter] << ", J = " << JTemp[iter] << "\n";
        // next iteration ======================================================
        iter++;
    }
    // summary screen output ===============================================
    std::cout << "\n\nSCVX Summary\nJ0=" << J[0] << "\n\n";
    for( int i=0; i<iter; i++ )
    {
        if( i%10 == 0 )
        {
            printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "Iteration","L","J","dL","dJ","rho","r","accept");
        }
        printf("%-12d%-12.6g%-12.6g%-12.3g%-12.3g%-12.3g%-12.3g%-12d\n",
               i+1,JTilde[i],JTemp[i],dL[i],dJ[i],rho[i],r[i],accept[i]);
    }
    return USucc;
}
