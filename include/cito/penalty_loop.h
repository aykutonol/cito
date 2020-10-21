/*! Penalty Loop */
/**
 *  \brief PenaltyLoop class consists of functions to run the penalty loop approach
 *
 *  This class defines functions that are usedthat are used to run the penalty loop
 *  approach that is used to generalize CITO. 
 *
 *  \author Aykut Onol
 */

#ifndef PENALTY_LOOP_H
#define PENALTY_LOOP_H

#include "scvx.h"

class PenaltyLoop
{
public:
    /// Constructor
    PenaltyLoop(const mjModel *m_, Params *cp_, SCVX *scvx_);
    /// Destructor
    ~PenaltyLoop();
    // This function executes the penalty loop algorithm and returns the optimal control trajectory
    eigMd solve(const eigMd &U0);

private:
    /// MuJoCo model
    const mjModel *m;
    /// Objects
    Params *cp;
    SCVX *scvx;
    /// Variables
    int iter = 0, maxIter = 2, lastAcceptedIter = 0;
    bool *accepted, *poseTolMet, *kMaxTolMet, *costTolMet, stop = false;
    double initPenalty = 0.1, reducePenaltyBy = 0.3, maxPenalty = 20.0,
           *penalty, *deltaPenalty,
           posTol = 0.3, rotTol = 1.0, kMaxTol = 0.1,
           posError0, rotError0, *posError, *rotError,
           *kAvg, *kMax;
    eigMd UPre, USol, UOpt;
    trajectory trajSol;
};

#endif //PENALTY_LOOP_H
