// ***** DESCRIPTION ***********************************************************
// PenaltyLoop class defines functions that are used to run the penalty loop
// approach that is used to generalize CITO.

#include "cito/penalty_loop.h"

// ***** CONSTRUCTOR ***********************************************************
PenaltyLoop::PenaltyLoop(const mjModel *m_, Params *cp_, Control *cc_, SCVX *scvx_) : m(m_), cp(cp_), cc(cc_), scvx(scvx_)
{
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
    // initialize penalty
    penalty[0] = initPenalty;
}
// ***** DESTRUCTOR ************************************************************
PenaltyLoop::~PenaltyLoop()
{
    delete[] accepted;
    delete[] poseTolMet;
    delete[] kMaxTolMet;
    delete[] costTolMet;
    delete[] penalty;
    delete[] deltaPenalty;
    delete[] posError;
    delete[] rotError;
    delete[] kAvg;
    delete[] kMax;
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
        printf("\033[3;33m\n\nPenalty iteration %d:\033[0m\n", iter);
        // set the weight on stiffness to the penalty value
        cp->weight[3] = penalty[iter];
        // refresh SCVX if not the first iteration
        if (iter > 0)
            scvx->refresh();
        // solve the SCVX problem
        USol = scvx->solveSCVX(UPre);
        // roll-out the dynamics and log the data
        trajSol = {};
        trajSol = scvx->runSimulation(USol, false, true, 1.);
        // get the initial pose error for normalization
        if (iter == 0)
            posError0 = (cp->desiredPos.head(2) - trajSol.X.col(0).segment(cp->controlJointDOF0, 2)).norm();
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
            printf("\t\033[0;33mPose tolerance met: pos. error: %f <= %f, rot. error: %f <= %f. Penalty change: %.3f.\033[0m\n",
                   posError[iter], posTol, rotError[iter], rotTol, deltaPenalty[iter]);
        }
        else
        {
            if (iter == 0)
                deltaPenalty[iter] = -initPenalty / 2.;
            else
                deltaPenalty[iter] = -fabs(deltaPenalty[iter - 1]) / 2.;
            printf("\t\033[0;33mPose tolerance not met: pos. error: %f > %f, rot. error: %f > %f. Penalty change: %.3f.\033[0m\n",
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
                deltaPenalty[iter] = 1.5;
            }
            printf("\t\033[0;33mMax. stiffness tolerance not met: %f > %f. Penalty change: %.3f.\033[0m\n",
                   kMax[iter], kMaxTol, deltaPenalty[iter]);
        }
        // apply the change
        penalty[iter + 1] = penalty[iter] + deltaPenalty[iter];
        // bound the penalty
        penalty[iter + 1] = std::max(penalty[iter + 1], 0.0);
        penalty[iter + 1] = std::min(penalty[iter + 1], maxPenalty);
        printf("\n\t\033[0;33mNext penalty value: %f.\033[0m\n\n", penalty[iter + 1]);
        // accept/reject the solution
        if (poseTolMet[iter] && (kAvg[iter] <= kAvg[lastAcceptedIter] || lastAcceptedIter == 0 || kAvg[iter] == 0))
        {
            accepted[iter] = true;
            lastAcceptedIter = iter;
            UPre = USol; // update the initial trajectory for next iteration
            UOpt = USol; // update the optimal solution
        }
        // post-process
        if (applyPP && accepted[iter] && kAvg[iter] <= kThresh)
        {
            eigMd UPost = postProcess(USol, trajSol.X);

            UPre = UPost;
            UOpt = UPost;
        }
        // check stopping condition
        if (poseTolMet[iter] && kMaxTolMet[iter])
            stop = true;
        else
        {
            iter++;
            std::cout << "\n\tPress any key to continue...\n";
            std::cin.ignore();
        }
    }
    return UOpt;
}

// postProcess: performs the post-process and returns the modified control trajectory
eigMd PenaltyLoop::postProcess(const eigMd &UPre, const eigMd &XPre)
{
    // apply the pulling controller
    eigMd UPost = pullControl(UPre, XPre);
    // apply the hill-climbing search
    // return the post-processed trajectory
    return UPost;
}

// pullControl: pulls virtually-active end effectors toward the corresponding contact geoms
eigMd PenaltyLoop::pullControl(const eigMd &UIn, const eigMd &XIn)
{
    eigMd UOut = UIn;
    // create & initialize MuJoCo data
    mjData *d = mj_makeData(m);
    mju_copy(d->qpos, m->key_qpos, m->nq);
    mj_forward(m, d);
    // intermediate variables
    eigVd pullVec(3), pullDir(3), pullFext(3), pullTau(cp->nu), sRbtPos(3), sEnvPos(3), nCS(3),
        qvelErr(cp->nv), qvelErrMass(cp->nv), dampTau(cp->nu);
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
            if (cp->taskType == 1)
                pullVec = sEnvPos + nCS * (pullControlShift / ((double)lastAcceptedIter + 1.)) - sRbtPos;
            else
                pullVec = sEnvPos - nCS * (pullControlShift / ((double)lastAcceptedIter + 1.)) - sRbtPos;
            pullDir = pullVec / pullVec.norm();
            // calculate the external pulling force acting on the end effector
            pullFext = pullControlKp * pullDir * UIn.col(i)(cp->nu + pair);
            std::cout << "\nt: " << d->time << " s\n\tpull dir: " << pullDir.transpose() << ", pull f: " << pullFext.transpose() << "\n\n";
            // get the translational Jacobian for the end-effector site
            mj_jacSite(m, d, JtSite.data(), NULL, cp->sites[pair][0]);
            JtRbt = JtSite.block(0, cp->dAct[0], 3, cp->nu);
            // project the pulling force onto the joint space
            pullTau += JtRbt.transpose() * pullFext;
        }
        // velocity tracking control
        qvelErr = XIn.col(i).tail(cp->nv) - XPull.col(i).tail(cp->nv);
        mj_mulM(m, d, qvelErrMass.data(), qvelErr.data());
        dampTau = dampControlKv * qvelErrMass.tail(cp->nu); // this assumes the actuated DOF are at the tail
        // apply the changes due to the pulling and velocity tracking controllers
        UOut.topRows(cp->nu).col(i) += pullTau + dampTau;
        // take a control step
        cc->takeStep(d, UOut.col(i), true, 1.);
    }
    // delete MuJoCo data
    mj_deleteData(d);
    // return the trajectory processed by the pulling controller
    return UOut;
}