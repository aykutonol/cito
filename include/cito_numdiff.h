/*! Numerical Differentiation */
/**
 *  \brief CitoNumDiff class consists of functions for numerical differentiation.
 *
 *  This class defines functions for numerical differentiation of the MuJoCo
 *  dynamics including the forces imposed by the contact model.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_NUMDIFF_H
#define CITO_NUMDIFF_H

#include "cito_control.h"

class CitoNumDiff
{
public:
    /// Constructor
    CitoNumDiff(const mjModel* model);
    /// Destructor
    ~CitoNumDiff() {}
    /// This function calculates derivatives of the state and control trajectories
    void linDyn(const mjData* dMain, const eigVd uMain, mjtNum* Fxd, mjtNum* Fud, double compensateBias);
    // This function performs fast finite-difference computation
    void worker(const mjData* dMain, mjtNum* deriv);

private:
    /// This function sets xNew to the integration of data given a control input
    void copyTakeStep(const mjData* dMain, const eigVd u, mjtNum* xNew, double compensateBias);
    /// This function calculates central differences by taking full steps
    void hardWorker(const mjData* dMain, const eigVd uMain, mjtNum* deriv, double compensateBias);
    /// MuJoCo model
    const mjModel* m;
    /// Perturbation for central differences
    double eps = 1e-8;
    /// Number of warm-up iterations for the worker function
    int nwarmup = 3;
    /// Variables for differentiation
    eigVm xNewTemp, xNewP, xNewN;
    eigVd uTemp;
    /// Objects
    CitoParams  cp;
    CitoControl cc;
};

#endif //CITO_NUMDIFF_H
