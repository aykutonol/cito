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
    int iter, maxIter = 1;
    bool *accepted, costTolMet = false, stop = false;
};

#endif //PENALTY_LOOP_H
