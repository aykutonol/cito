// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_NUMDIFF class consists of functions for numerical differentiation of the
// MuJoCo dynamics including the contact forces imposed by the contact model.

#include "cito_control.h"

#ifndef CITO_NUMDIFF_H
#define CITO_NUMDIFF_H

class CitoNumDiff
{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
    CitoNumDiff(const mjModel* model);
    ~CitoNumDiff() {}
    // ***** FUNCTIONS *************************************************************
    void linDyn(const mjData* dmain, const ctrlVec_t umain, mjtNum* Fxd, mjtNum* Fud);

private:
    // ***** FUNCTIONS *************************************************************
    void copyTakeStep(const mjData* dmain, const ctrlVec_t u, mjtNum* newXd);
    void hardWorker(const mjData* dmain, const ctrlVec_t umain, mjtNum* deriv);
    // ***** PARAMETERS ************************************************************
    const mjModel* m;
    // perturbation
    double eps = 1e-6;
    // variables for differentiation
    stateVec_t newXtemp, newXp, newXn;
    ctrlVec_t  utemp;
    // ***** OBJECTS ***************************************************************
    CitoControl cc;
};

#endif //CITO_NUMDIFF_H