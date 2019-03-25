#include "cito_scvx.h"

// ***** DESCRIPTION ***********************************************************
// CitoSCvx class defines functions that are used to roll-out the dynamics,
// evaluate the cost, and execute the SCvx algorithm.

// ***** CONSTRUCTOR ***********************************************************
CitoSCvx::CitoSCvx(const mjModel* model) : m(model), cp(model), cc(model), nd(model), sq(model)
{
    // initialize Eigen variables
    desiredPos.resize(6);       finalPos.resize(6);
    desiredVel.resize(m->nv);   finalVel.resize(m->nv);
    KCon.resize(cp.nPair,cp.N);
    // waypoint parameters
    desiredPosWP.resize(6);     wayPtPos.resize(6);
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
    XSucc.resize(cp.n,cp.N+1);  dX.resize(cp.n,cp.N+1);     XTilde.resize(cp.n,cp.N+1);
    USucc.resize(cp.m,cp.N);    UTemp.resize(cp.m,cp.N);    dU.resize(cp.m,cp.N);
    Fx.resize(cp.n*cp.n,cp.N);  Fu.resize(cp.n*cp.m,cp.N);
    // read cost function weights
    weight[0] = params["w1"].as<double>();
    weight[1] = params["w2"].as<double>();
    weight[2] = params["w3"].as<double>();
    weight[3] = params["w4"].as<double>();
}

// ***** FUNCTIONS *************************************************************
// getCost: returns the nonlinear cost given state and control trajectories
double CitoSCvx::getCost(const eigMm X, const eigMd U)
{
    // final cost
    finalPos = X.col(cp.N).segment(controlJointDOF0,6);
    wayPtPos = X.col(tWP).segment(controlJointDOF0,6);
    finalVel = X.col(cp.N).tail(m->nv);
    Jf = 0.5*(weight[0]*(desiredPos.head(2)-finalPos.head(2)).squaredNorm()+
              weight[1]*(desiredPos.tail(4)-finalPos.tail(4)).squaredNorm()+
              weight[0]*(desiredPosWP.head(2)-wayPtPos.head(2)).squaredNorm()+
              weight[1]*(desiredPosWP.tail(4)-wayPtPos.tail(4)).squaredNorm()+
              weight[2]*(desiredVel - finalVel).squaredNorm());
    // integrated cost
    KConSN = 0;
    for( int i=0; i<cp.N; i++ )
    {
        KCon.col(i).setZero();
        for( int j=0; j<cp.nPair; j++ )
        {
            KCon.col(i)(j) = U.col(i)[m->nu+j];
        }
        KConSN += KCon.col(i).squaredNorm();
    }
    Ji = 0.5*weight[3]*KConSN;
    // total cost
    Jt = Jf + Ji;
    return Jt;
}

// runSimulation: rolls-out and linearizes the dynamics given control trajectory
trajectory CitoSCvx::runSimulation(const eigMd U, bool linearize, bool save)
{
    // make mjData
    mjData* d = NULL;
    d = mj_makeData(m);
    // initialize d
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    cc.setControl(d, U.col(0));
    // rollout (and linearize) the dynamics
    for( int i=0; i<cp.N; i++ )
    {
        // get the current state values
        XSucc.col(i).setZero();
        XSucc.col(i) = cc.getState(d);
        // linearization
        if( linearize )
        {
            Fx.col(i).setZero(); Fu.col(i).setZero();
            nd.linDyn(d, U.col(i), Fx.col(i).data(), Fu.col(i).data());
        }
        // take tc/dt steps
        cc.takeStep(d, U.col(i), save);
    }
    XSucc.col(cp.N).setZero();
    XSucc.col(cp.N) = cc.getState(d);
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
eigMd CitoSCvx::solveSCvx(const eigMd U0)
{
    // initialize USucc for the first succession
    for( int i=0; i<cp.N; i++ ) { USucc.col(i) = U0.col(i); }
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
        if( iter == 0 ) { J[iter] = this->getCost(trajS.X, USucc); }
        // convex optimization =================================================
        double *dTraj = new double[cp.nTraj];
        std::cout << "INFO: QP solver in progress\n\n";
        auto tQPStart = std::chrono::system_clock::now();
        sq.solveCvx(dTraj, r[iter], trajS.X, USucc, trajS.Fx, trajS.Fu, cc.isJFree, cc.isAFree,
                    cc.qposLB, cc.qposUB, cc.tauLB, cc.tauUB);
        auto tQPEnd = std::chrono::system_clock::now();
        std::cout << "\nINFO: QP solver took " << std::chrono::duration<double>(tQPEnd-tQPStart).count() << " s \n\n";
        // apply the change
        for( int i=0; i<cp.N+1; i++ )
        {
            // states
            dX.col(i).setZero(); XTilde.col(i).setZero();
            for( int j=0; j<cp.n; j++ )
            {
                dX.col(i)[j] = dTraj[i*cp.n+j];
            }
            XTilde.col(i) = trajS.X.col(i) + dX.col(i);
            // controls
            if( i < cp.N )
            {
                dU.col(i).setZero(); UTemp.col(i).setZero();
                for( int j=0; j<cp.m; j++ )
                {
                    dU.col(i)[j] = dTraj[(cp.N+1)*cp.n+i*cp.m+j];
                }
                UTemp.col(i) = USucc.col(i) + dU.col(i);
            }
        }
        // evaluate the dynamics for the change and get the cost values ========
        trajTemp = {};
        trajTemp = this->runSimulation(UTemp, false, false);
        // get the linear and nonlinear costs
        JTilde[iter] = this->getCost(XTilde, UTemp);
        JTemp[iter]  = this->getCost(trajTemp.X, UTemp);
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
            for( int i=0; i<cp.N; i++ )
            {
                USucc.col(i) = UTemp.col(i);
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
            std::cout << "\n\n\tINFO: Maximum number of iterations reached.\n\n";
        }
        if( dLTolMet )
        {
            stop = true;
            std::cout << "\n\n\tINFO: |dL| = |" << dL[iter] << "| < dLTol = " << dLTol << "\n\n";
        }
        // screen output for the iteration =====================================
        std::cout << "Actual:\nFinal pos: " << trajTemp.X.col(cp.N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << trajTemp.X.col(cp.N).tail(m->nv).transpose() << "\n";
        std::cout << "Predicted:\nFinal pos: " << XTilde.col(cp.N).head(m->nv).transpose() << "\n";
        std::cout << "Final vel: " << XTilde.col(cp.N).tail(m->nv).transpose() << "\n";
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
