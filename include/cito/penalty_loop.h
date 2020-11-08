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
    PenaltyLoop(const mjModel *m_, Params *cp_, Control *cc_, SCVX *scvx_);
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
    Control *cc;
    /// Variables
    int iter = 0, maxIter = 10, lastAcceptedIter = 0,
        iterHCS, maxIterHCS = 100;
    bool stop = false, *accepted, *poseTolMet, *kMaxTolMet, *costTolMet,
         applyPP = true, *acceptPP;
    double initPenalty = 0.1, reducePenaltyBy = 0.3, maxPenalty = 20.0,
           *penalty, *deltaPenalty,
           posTol = 0.3, rotTol = 1., kMaxTol = 0.1, kThresh = 2.5,
           posError0, rotError0, *posError, *rotError,
           *kAvg, *kMax,
           pullControlShift = 0.05, pullControlKp = 5., dampControlKv = 1.,
           *costHCS, alphaHCS = 10., dCostHCS, lastCostHCS,
           kAvgBeforeHCS, kMaxBeforeHCS, *kAvgReduction, *kMaxReduction;
    eigMd UPre, USol, UOpt;
    trajectory trajSol, trajHCS;
    // This function performs the post-process and returns the modified control trajectory
    eigMd postProcess(const eigMd &UPre, const eigMd &XPre);
    // This function pulls virtually-active end effectors toward the corresponding contact geoms
    trajectory pullControl(const eigMd &UIn, const eigMd &XIn);
    // This function reduces the excessive virtual stiffness values after applying the pull control
    eigMd hillClimbSearch(eigMd &UPulled, const eigMd &XPulled);
};

#endif //PENALTY_LOOP_H
