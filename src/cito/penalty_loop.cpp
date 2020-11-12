// ***** DESCRIPTION ***********************************************************
// PenaltyLoop class defines functions that are used to run the penalty loop
// approach that is used to generalize CITO.

#include "cito/penalty_loop.h"

// ***** CONSTRUCTOR ***********************************************************
PenaltyLoop::PenaltyLoop(const mjModel *m_, Params *cp_, Control *cc_, SCVX *scvx_) : m(m_), cp(cp_), cc(cc_), scvx(scvx_)
{
    // get penalty loop parameters
    YAML::Node paramPenalty = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/penalty_loop.yaml");
    maxIter = paramPenalty["maxIter"].as<int>();
    maxPenalty = paramPenalty["maxPenalty"].as<double>();
    initPenalty = paramPenalty["initPenalty"].as<double>();
    incPenaltyFixed = paramPenalty["incPenaltyFixed"].as<double>();
    decPenaltyRatio = paramPenalty["decPenaltyRatio"].as<double>();
    posTol = paramPenalty["posTol"].as<double>();
    rotTol = paramPenalty["rotTol"].as<double>();
    kMaxTol = paramPenalty["kMaxTol"].as<double>();
    kThresh = paramPenalty["kThresh"].as<double>();
    pullControlKp = paramPenalty["Kp"].as<double>();
    dampControlKv = paramPenalty["Kv"].as<double>();
    pullControlShift = paramPenalty["shift"].as<double>();
    alphaHCS = paramPenalty["alpha"].as<double>();
    maxIterHCS = paramPenalty["maxIterHCS"].as<int>();
    // create dynamic arrays
    accepted = new bool[maxIter]();   // initialize with false
    acceptPP = new bool[maxIter]();   // initialize with false
    poseTolMet = new bool[maxIter](); // initialize with false
    kMaxTolMet = new bool[maxIter](); // initialize with false
    costTolMet = new bool[maxIter](); // initialize with false
    penalty = new double[maxIter + 1];
    deltaPenalty = new double[maxIter];
    posError = new double[maxIter];
    rotError = new double[maxIter];
    kAvg = new double[maxIter];
    kMax = new double[maxIter];
    kAvgReduction = new double[maxIter];
    kMaxReduction = new double[maxIter];
    costHCS = new double[maxIterHCS];
    // initialize penalty
    penalty[0] = initPenalty;
}
// ***** DESTRUCTOR ************************************************************
PenaltyLoop::~PenaltyLoop()
{
    delete[] accepted;
    delete[] acceptPP;
    delete[] poseTolMet;
    delete[] kMaxTolMet;
    delete[] costTolMet;
    delete[] penalty;
    delete[] deltaPenalty;
    delete[] posError;
    delete[] rotError;
    delete[] kAvg;
    delete[] kMax;
    delete[] kAvgReduction;
    delete[] kMaxReduction;
    delete[] costHCS;
}

// ***** FUNCTIONS *************************************************************
// Solve: executes the penalty loop algorithm and returns the optimal control trajectory
eigMd PenaltyLoop::solve(const eigMd &U0)
{
    // copy the initial control trajectory
    UPre = U0;
    UOpt = U0;
    // run the penalty loop algorithm
    while (iter < maxIter && !stop)
    {
        // penalty loop info
        printf("\033[3;33m\n\nPenalty iteration %d:\033[0m\n", iter + 1);
        // set the weight on stiffness to the penalty value
        cp->weight[3] = penalty[iter];
        // refresh SCVX if not the first iteration
        if (iter > 0)
            scvx->refresh();
        // solve the SCVX problem
        USol = scvx->solveSCVX(UPre);
        // roll-out the dynamics and log the data
        trajSol = {};
        trajSol = scvx->runSimulation(USol, false, 1, 1.);
        // get the initial pose error for normalization
        if (iter == 0)
            posError0 = std::max(1e-2, (cp->desiredPos.head(2) - trajSol.X.col(0).segment(cp->controlJointDOF0, 2)).norm());
        // get average and maximum stiffness
        kAvg[iter] = USol.bottomRows(cp->nPair).sum() / (cp->nPair * cp->N);
        kMax[iter] = USol.bottomRows(cp->nPair).maxCoeff();
        // get the position (normalized) and rotation errors
        posError[iter] = (cp->desiredPos.head(2) - trajSol.X.col(cp->N).segment(cp->controlJointDOF0, 2)).norm() / posError0;
        rotError[iter] = (cp->desiredPos.tail(3) - trajSol.X.col(cp->N).segment(cp->controlJointDOF0 + 3, 3)).norm();
        // check if the pose tolerance is satisfied
        if (posError[iter] <= posTol && rotError[iter] <= rotTol)
        {
            poseTolMet[iter] = true;
            deltaPenalty[iter] = 0.;
            printf("\t\033[0;33mPose tolerance met: pos. error: %f <= %f and rot. error: %f <= %f. Penalty change: %.3f.\033[0m\n",
                   posError[iter], posTol, rotError[iter], rotTol, deltaPenalty[iter]);
        }
        else
        {
            if (iter == 0)
                deltaPenalty[iter] = -initPenalty / 2.;
            else
                deltaPenalty[iter] = -fabs(deltaPenalty[iter - 1]) * decPenaltyRatio;
            printf("\n\t\033[0;33mPose tolerance not met: pos. error: %f > %f or rot. error: %f > %f. Penalty change: %.3f.\033[0m\n",
                   posError[iter], posTol, rotError[iter], rotTol, deltaPenalty[iter]);
        }
        // check if the maximum stiffness tolerance is satisfied
        if (kMax[iter] <= kMaxTol)
        {
            kMaxTolMet[iter] = true;
            printf("\t\033[0;33mMax. stiffness tolerance met: %f <= %f. Penalty change: %.3f.\033[0m\n",
                   kMax[iter], kMaxTol, deltaPenalty[iter]);
        }
        else
        {
            if (poseTolMet[iter] || (!poseTolMet[iter] && penalty[iter] == 0.0))
            {
                deltaPenalty[iter] = incPenaltyFixed;
            }
            printf("\t\033[0;33mMax. stiffness tolerance not met: %f > %f. Penalty change: %.3f.\033[0m\n",
                   kMax[iter], kMaxTol, deltaPenalty[iter]);
        }
        // apply the change
        penalty[iter + 1] = penalty[iter] + deltaPenalty[iter];
        // bound the penalty
        penalty[iter + 1] = std::max(penalty[iter + 1], 0.0);
        penalty[iter + 1] = std::min(penalty[iter + 1], maxPenalty);
        printf("\t\033[0;33mNext penalty value: %f.\033[0m\n\n", penalty[iter + 1]);
        // accept/reject the solution
        if (poseTolMet[iter] && (kAvg[iter] <= kAvg[lastAcceptedIter] || lastAcceptedIter == 0 || kAvg[iter] == 0))
        {
            printf("\t\033[0;32mPenalty iteration %d accepted.\033[0m\n", iter + 1);
            accepted[iter] = true;
            lastAcceptedIter = iter;
            UPre = USol; // update the initial trajectory for next iteration
            UOpt = USol; // update the optimal solution
        }
        else
            printf("\t\033[0;31mPenalty iteration %d rejected.\033[0m\n", iter + 1);
        // post-process
        if (applyPP && accepted[iter] && kAvg[iter] <= kThresh)
        {
            printf("\t\033[0;33mPerforming the post-process...\033[0m\n");
            eigMd UPost = postProcess(USol, trajSol.X);
            trajectory trajPost = scvx->runSimulation(UPost, false, 1, 1.);
            UPre = UPost;
            UOpt = UPost;
        }
        // check stopping condition
        if (poseTolMet[iter] && kMaxTolMet[iter])
            stop = true;
        else
        {
            iter++;
            std::cout << "\n\tPenalty iteration " << iter << " completed. Press any key to continue...\n";
            std::cin.ignore();
        }
    }
    return UOpt;
}

// postProcess: performs the post-process and returns the modified control trajectory
eigMd PenaltyLoop::postProcess(const eigMd &UPre, const eigMd &XPre)
{
    // apply the pulling controller
    trajectory trajPull = pullControl(UPre, XPre);
    // apply the hill-climbing search
    eigMd UPost = hillClimbSearch(trajPull.U, trajPull.X);
    // return the post-processed trajectory
    return UPost;
}

// pullControl: pulls virtually-active end effectors toward the corresponding contact geoms
trajectory PenaltyLoop::pullControl(const eigMd &UIn, const eigMd &XIn)
{
    trajectory trajOut;
    trajOut.U = UIn;
    // create & initialize MuJoCo data
    mjData *d = mj_makeData(m);
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    // intermediate variables
    eigVd pullVec(3), pullDir(3), pullFext(3), pullTau(cp->nu), sRbtPos(3), sEnvPos(3), nCS(3),
        qvelError(cp->nv), qvelErrorMass(cp->nv), dampTau(cp->nu);
    eigMd XPull(cp->n, cp->N + 1);
    XPull.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> JtSite(3, cp->nv), JtRbt(3, cp->nu);
    // pass the trajectory through the pulling controller
    for (int i = 0; i < cp->N; i++)
    {
        // get the current state
        XPull.col(i) = cp->getState(d);
        // set the pulling joint force to zero
        pullTau.setZero();
        // induce the pulling behavior
        for (int pair = 0; pair < cp->nPair; pair++)
        {
            // get the site positions on the robot and in the environment
            mju_copy3(sRbtPos.data(), d->site_xpos + 3 * cp->sites[pair][0]);
            mju_copy3(sEnvPos.data(), d->site_xpos + 3 * cp->sites[pair][1]);
            // update contact surface normal
            mju_rotVecMat(nCS.data(), cp->unit_x, d->site_xmat + 9 * cp->sites[pair][1]);
            // shift the pulling target out along the surface normal in the environment
            // normal pointing outwards for locomotion and inwards for manipulation tasks
            reducedShift = pullControlShift / pow(((double)lastAcceptedIter + 1.), 1.);
            if (cp->taskType == 1)
                pullVec = sEnvPos - sRbtPos + nCS * reducedShift;
            else
                pullVec = sEnvPos - sRbtPos - nCS * reducedShift;
            pullDir = pullVec / pullVec.norm();
            // calculate the external pulling force acting on the end effector
            pullFext = pullControlKp * pullDir * UIn.col(i)(cp->nu + pair);
            if (UIn.col(i)(cp->nu + pair) > 1e-3)
                std::cout << "\nt: " << d->time << " s, pair: " << pair
                          << "\n\tpull dir: " << pullDir.transpose() << ", pull f: " << pullFext.transpose()
                          << "\n\tsite rbt: " << sRbtPos.transpose() << ", site env: " << sEnvPos.transpose() << "\n\n";
            // get the translational Jacobian for the end-effector site
            mj_jacSite(m, d, JtSite.data(), NULL, cp->sites[pair][0]);
            JtRbt = JtSite.block(0, cp->dAct[0], 3, cp->nu);
            // project the pulling force onto the joint space
            pullTau += JtRbt.transpose() * pullFext;
        }
        // velocity tracking control
        qvelError = XIn.col(i).tail(cp->nv) - XPull.col(i).tail(cp->nv);
        mj_mulM(m, d, qvelErrorMass.data(), qvelError.data());
        dampTau = dampControlKv * qvelErrorMass.tail(cp->nu); // this assumes the actuated DOF are at the tail
        // apply the changes due to the pulling and velocity tracking controllers
        trajOut.U.topRows(cp->nu).col(i) += pullTau + dampTau;
        // take a control step
        cc->takeStep(d, trajOut.U.col(i), 0, 1.);
    }
    // get the final state
    XPull.col(cp->N) = cp->getState(d);
    trajOut.X = XPull;
    // delete MuJoCo data
    mj_deleteData(d);
    // return the trajectory processed by the pulling controller
    return trajOut;
}

// hillClimbSearch reduces the excessive virtual stiffness values after applying the pull control
eigMd PenaltyLoop::hillClimbSearch(eigMd &UPulled, const eigMd &XPulled)
{
    eigMd UHCS = UPulled;
    // get the initial cost, avg. and max. stiffness values
    costHCS[0] = (cp->desiredPos.head(3) - XPulled.col(cp->N).segment(cp->controlJointDOF0, 3)).norm() +
                 (cp->weight[1] / cp->weight[0]) * (cp->desiredPos.tail(3) - XPulled.col(cp->N).segment(cp->controlJointDOF0 + 3, 3)).norm();
    lastCostHCS = costHCS[0];
    kAvgBeforeHCS = UPulled.bottomRows(cp->nPair).sum() / (cp->nPair * cp->N);
    kMaxBeforeHCS = UPulled.bottomRows(cp->nPair).maxCoeff();
    // run the hill-climbing search
    iterHCS = 0;
    std::cout << "Initial cost: " << costHCS[0]
              << ", stiffness trajectory:\n"
              << UPulled.bottomRows(cp->nPair) << "\n";
    printf("\033[0;33mStarting hill-climbing search\033[0m\n");
    for (int pair = 0; pair < cp->nPair; pair++)
    {
        for (int i = 0; i < cp->N; i++)
        {
            if (UPulled(cp->nu + pair, i) > 0.)
            {
                dCostHCS = 1e-2;
                UHCS = UPulled;
                while (dCostHCS > 1e-4 && iterHCS < maxIterHCS && costHCS[iterHCS] > 1e-2)
                {
                    iterHCS++;
                    // UHCS(cp->nu + pair, i) -= alphaHCS * dCostHCS;
                    UHCS(cp->nu + pair, i) -= 0.1;
                    UHCS(cp->nu + pair, i) = std::max(0., UHCS(cp->nu + pair, i));
                    std::cout << "\tApplied change: " << -alphaHCS * dCostHCS
                              << " on pair " << pair
                              << " at time step " << i
                              << "\nStiffness values:\n"
                              << UHCS.bottomRows(cp->nPair) << "\n";
                    trajHCS = scvx->runSimulation(UHCS, false, 0, 1.);
                    costHCS[iterHCS] = (cp->desiredPos.head(3) - trajHCS.X.col(cp->N).segment(cp->controlJointDOF0, 3)).norm() +
                                       (cp->weight[1] / cp->weight[0]) * (cp->desiredPos.tail(3) - trajHCS.X.col(cp->N).segment(cp->controlJointDOF0 + 3, 3)).norm();
                    dCostHCS = lastCostHCS - costHCS[iterHCS];
                    if (dCostHCS > 0.)
                    {
                        UPulled = UHCS;
                        lastCostHCS = costHCS[iterHCS];
                    }
                    printf("\tHCS iter: %d, cost: %f, dcost: %f, alpha: %f\n\n",
                           iterHCS, costHCS[iterHCS], dCostHCS, alphaHCS);
                }
            }
        }
    }
    if (iterHCS == 0)
        printf("\033[0;33mInitial cost is smaller than 0.01.\033[0m\n");
    // calculate the reduction of avg. and max. stiffness
    kAvgReduction[iter] = kAvgBeforeHCS - UPulled.bottomRows(cp->nPair).sum() / (cp->nPair * cp->N);
    kMaxReduction[iter] = kMaxBeforeHCS - UPulled.bottomRows(cp->nPair).maxCoeff();
    // return the output of the hill-climbing search
    return UPulled;
}