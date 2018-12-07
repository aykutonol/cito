// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_NUMDIFF class consists of functions for numerical differentiation of the
// MuJoCo dynamics including the contact forces imposed by the contact model.

#include "cito_params.h"
#include "cito_control.h"

#ifndef CITO_NUMDIFF_H
#define CITO_NUMDIFF_H

class CitoNumdiff
{
private:
    const mjModel*  m;
    const mjData*   dmain;
    const ctrlVec_t umain;
    // perturbation
    double eps = 1e-6;
    // variables for differentiation
    stateVec_t newXtemp, newXp, newXn;
    ctrlVec_t up;
    // custom objects
    CitoControl cc;
    void CitoNumDiff::copyTakeStep(mjtNum* newX, ctrlVec_t u);
    void CitoNumDiff::hardWorker(mjtNum* deriv);

public:
    CitoNumDiff(const mjModel* model, const mjData* data, const ctrlVec_t control);
    ~CitoNumDiff() {}
    stateVec_t CitoNumDiff::getState(const mjData* d);
    void CitoNumDiff::linDyn(mjtNum* Fxd, mjtNum* Fud);
}

#endif //CITO_NUMDIFF_H
