// ***** DESCRIPTION ***********************************************************
// PenaltyLoop class defines functions that are used to run the penalty loop
// approach that is used to generalize CITO. 

#include "cito/penalty_loop.h"

// ***** CONSTRUCTOR ***********************************************************
PenaltyLoop::PenaltyLoop(const mjModel *m_, Params *cp_, SCVX *scvx_) : m(m_), cp(cp_), scvx(scvx_)
{
    accepted = new bool[maxIter];
}
// ***** DESTRUCTOR ************************************************************
PenaltyLoop::~PenaltyLoop()
{
    delete[] accepted;
}

// ***** FUNCTIONS *************************************************************
// Solve: executes the penalty loop algorithm and returns the optimal control trajectory
eigMd PenaltyLoop::solve(const eigMd &U0)
{
    eigMd Uopt;
    for (int i=0; i<maxIter; i++)
    {
        Uopt = scvx->solveSCVX(U0);
    }
    return Uopt;
}