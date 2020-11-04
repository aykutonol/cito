// ***** DESCRIPTION ***********************************************************
// PenaltyLoop class defines functions that are used to run the penalty loop
// approach that is used to generalize CITO.

#include "cito/penalty_loop.h"

// ***** CONSTRUCTOR ***********************************************************
PenaltyLoop::PenaltyLoop(const mjModel *m_, Params *cp_, SCVX *scvx_) : m(m_), cp(cp_), scvx(scvx_)
{
    accepted = new bool[maxIter]();   // initialize with false
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
        // roll-out the dynamics
        trajSol = {};
        trajSol = scvx->runSimulation(USol, false, false, 1);
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