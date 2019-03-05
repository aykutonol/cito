/*! Numerical Differentiation */
/**
 *  \brief CitoSCvx class consists of functions for numerical differentiation.
 *
 *  This class defines functions for numerical differentiation of the MuJoCo
 *  dynamics including the forces imposed by the contact model.
 *
 *  \author Aykut Onol
 */

#include "cito_control.h"

#ifndef CITO_NUMDIFF_H
#define CITO_NUMDIFF_H

class CitoNumDiff
{
public:
    /// Constructor
    CitoNumDiff(const mjModel* model);
    /// Destructor
    ~CitoNumDiff() {}
    /// This function calculates derivatives of the state and control trajectories
    void linDyn(const mjData* dmain, const ctrlVec_t umain, mjtNum* Fxd, mjtNum* Fud);

private:
    /// This function sets newXd to the integration of data given a control input
    void copyTakeStep(const mjData* dmain, const ctrlVec_t u, mjtNum* newX);
    /// This function calculates central differences by taking full steps
    void hardWorker(const mjData* dmain, const ctrlVec_t umain, mjtNum* deriv);
    /// MuJoCo model
    const mjModel* m;
    /// Perturbation for central differences
    double eps = 1e-6;
    /// Variables for differentiation
    stateVec_t newXtemp, newXp, newXn;
    ctrlVec_t  utemp;
    /// Control object
    CitoControl cc;
};

#endif //CITO_NUMDIFF_H
