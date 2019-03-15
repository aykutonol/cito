#include "cito_scvx.h"

// ***** DESCRIPTION ***********************************************************
// CitoSCvx class defines functions that are used to roll-out the dynamics,
// evaluate the cost, and execute the SCvx algorithm.

// ***** CONSTRUCTOR ***********************************************************
CitoSCvx::CitoSCvx(const mjModel* model) : m(model), cc(model), nd(model)
{
    // read task parameters
    YAML::Node paramTask = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/task.yaml");
    std::vector<double> desiredPosInput = { paramTask["desiredFinalPos"].as<std::vector<double>>() };
    std::vector<double> desiredVelInput = { paramTask["desiredFinalVel"].as<std::vector<double>>() };
    desiredPos = Eigen::Map<Eigen::Matrix<double,  6, 1>>(desiredPosInput.data(), desiredPosInput.size());
    desiredVel = Eigen::Map<Eigen::Matrix<double, NV, 1>>(desiredVelInput.data(), desiredVelInput.size());
    controlJointDOF0 = paramTask["controlJointDOF0"].as<int>();
    // read SCvx parameters
    YAML::Node paramSCvx = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/scvx.yaml");
    // create new arrays for the max. number of iterations
    beta_expand = paramSCvx["beta_expand"].as<double>();
    beta_shrink = paramSCvx["beta_shrink"].as<double>();
    maxIter = paramSCvx["maxIter"].as<int>();
    dLTol = paramSCvx["dLTol"].as<double>();
    rho0 = paramSCvx["rho0"].as<double>();
    rho1 = paramSCvx["rho1"].as<double>();
    rho2 = paramSCvx["rho2"].as<double>();
    rMin = paramSCvx["rMin"].as<double>();
    rMax = paramSCvx["rMax"].as<double>();
    J      = new double[maxIter+1];
    JTemp  = new double[maxIter+1];
    JTilde = new double[maxIter+1];
    dJ     = new double[maxIter+1];
    dL     = new double[maxIter+1];
    rho    = new double[maxIter+1];
    r      = new double[maxIter+1];
    accept = new bool[maxIter];
    r[0]  = paramSCvx["r0"].as<double>();        // initial trust-region radius
    // set bounds
    cc.getBounds();
    // resize trajectories
    XSucc.resize(NTS+1);    dX.resize(NTS+1);   XTilde.resize(NTS+1);
    USucc.resize(NTS);      UTemp.resize(NTS);  dU.resize(NTS);
    Fx.resize(NTS);         Fu.resize(NTS);
    KCon.resize(NTS);
    // read cost function weights
    weight[0] = paramTask["w1"].as<double>();
    weight[1] = paramTask["w2"].as<double>();
    weight[2] = paramTask["w3"].as<double>();
    weight[3] = paramTask["w4"].as<double>();
}

// ***** FUNCTIONS *************************************************************
// getCost: returns the nonlinear cost given control trajectory and final state
double CitoSCvx::getCost(stateVec XFinal, const ctrlTraj U)
{
    // terminal cost
    for( int i=0; i<6; i++ )
    {
        finalPos[i] = XFinal[controlJointDOF0 + i];
    }
    for( int i=0; i<NV; i++ )
    {
        finalVel[i] = XFinal[controlJointDOF0 + NV + i];
    }
    Jf = 0.5*(weight[0]*(desiredPos.block<2,1>(0,0)-finalPos.block<2,1>(0,0)).squaredNorm()+
              weight[1]*(desiredPos.block<4,1>(2,0)-finalPos.block<4,1>(2,0)).squaredNorm()+
              weight[2]*(desiredVel - finalVel).squaredNorm());
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
    Ji = 0.5*weight[3]*KConSN;
    // total cost
    Jt = Jf + Ji;
    return Jt;
}

// runSimulation: rolls-out and linearizes the dynamics given control trajectory
trajectory CitoSCvx::runSimulation(const ctrlTraj U, bool linearize, bool save)
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
ctrlTraj CitoSCvx::solveSCvx(const ctrlTraj U0)
{
    // initialize USucc for the first succession
    for( int i=0; i<NTS; i++ ) { USucc[i] = U0[i]; }
    // start the SCvx algorithm
    int iter = 0;
    while( !stop )
    {
        std::cout << "Iteration " << iter+1 << ":" << '\n';
        // simulation and convexification ======================================
        if( iter == 0 || accept[iter-1] )
        {
            std::cout << "INFO: convexification in progress\n";
            auto tDiffStart = std::chrono::system_clock::now();
            trajS = {};
            trajS = this->runSimulation(USucc, true, false);
            auto tDiffEnd = std::chrono::system_clock::now();
            std::cout << "INFO: convexification took " << std::chrono::duration<double>(tDiffEnd-tDiffStart).count() << " s \n";
        }
        // get the nonlinear cost if the first iteration
        if( iter == 0 ) { J[iter] = this->getCost(trajS.X[NTS], USucc); }
        // convex optimization =================================================
        double *dTraj = new double[NTRAJ];
        std::cout << "INFO: QP solver in progress\n\n";
        auto tQPStart = std::chrono::system_clock::now();
        sq.solveCvx(dTraj, r[iter], trajS.X, USucc, trajS.Fx, trajS.Fu, cc.isJFree, cc.isAFree,
                    cc.qposLB, cc.qposUB, cc.tauLB, cc.tauUB);
        auto tQPEnd = std::chrono::system_clock::now();
        std::cout << "\nINFO: QP solver took " << std::chrono::duration<double>(tQPEnd-tQPStart).count() << " s \n\n";
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
        JTilde[iter] = this->getCost(XTilde[NTS], UTemp);
        JTemp[iter]  = this->getCost(trajTemp.X[NTS], UTemp);
        // similarity measure ==================================================
        dJ[iter] = J[iter] - JTemp[iter];
        dL[iter] = J[iter] - JTilde[iter];
        rho[iter] = dJ[iter]/dL[iter];
        if( fabs(dL[iter]) < dLTol )
        {
            dLTolMet = 1;
        }
        // accept or reject the solution =======================================
        // reject
        if( rho[iter]<=rho0 || (dL[iter]<0 && dJ[iter]<0) )
        {
            accept[iter] = false;
            r[iter+1] = r[iter]/beta_shrink;
            J[iter+1] = J[iter];
        }
        else { accept[iter] = true; }
        // accept
        if( accept[iter] )
        {
            J[iter+1] = JTemp[iter];
            for( int i=0; i<NTS; i++ )
            {
                USucc[i] = UTemp[i];
            }
            if( rho[iter] < rho1 )
            { r[iter+1] = r[iter]/beta_shrink;  }
            else if( rho[iter]>=rho1 && rho[iter]<rho2 )
            { r[iter+1] = r[iter];        }
            else if( rho[iter]>=rho2 )
            { r[iter+1] = r[iter]*beta_expand;   }
        }
        // bound the trust region radius r
        r[iter+1] = std::max(r[iter+1], rMin);
        r[iter+1] = std::min(r[iter+1], rMax);
        // stopping criteria check =============================================
        if( iter+1 == maxIter )
        {
            stop = true;
            std::cout << "\n\n\tINFO: Maximum number of iterations reached\n\n";
        }
        if( dLTolMet )
        {
            stop = true;
            std::cout << "\n\n\tINFO: dL = " << fabs(dL[iter]) << " < dLTol = " << dLTol << "\n\n";
        }
        // screen output for the iteration =====================================
        std::cout << "Actual:\nFinal pos: " << trajTemp.X[NTS].block<NV,1>(0,0).transpose() << "\n";
        std::cout << "Final vel: " << trajTemp.X[NTS].block<NV,1>(NV,0).transpose() << "\n";
        std::cout << "Predicted:\nFinal pos: " << XTilde[NTS].block<NV,1>(0,0).transpose() << "\n";
        std::cout << "Final vel: " << XTilde[NTS].block<NV,1>(NV,0).transpose() << "\n";
        std::cout << "J = " << JTemp[iter] << ", JTilde = " << JTilde[iter] << "\n\n\n";
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
