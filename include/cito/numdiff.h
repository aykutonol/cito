/*! Numerical Differentiation */
/**
 *  \brief NumDiff consists of methods for numerical differentiation.
 *
 *  This class defines functions for numerical differentiation of the MuJoCo
 *  dynamics including the forces imposed by the contact model.
 *
 *  \author Aykut Onol
 */

#ifndef NUMDIFF_H
#define NUMDIFF_H

#include "cito/control.h"

class NumDiff
{
public:
    /// Constructor
    NumDiff(const mjModel *m_, Params *cp_, Control *cc_);
    /// Destructor
    ~NumDiff() {}
    /// This function calculates derivatives of the state and control trajectories
    void linDyn(const mjData *dMain, const eigVd &uMain, double *Fxd, double *Fud, double compensateBias);

private:
    /// This function sets xNew to the integration of data given a control input
    void copyTakeStep(const mjData *dMain, const eigVd &u, double *xNew, double compensateBias);
    /// This function calculates central differences by taking full steps
    void hardWorker(const mjData *dMain, const eigVd &uMain, double *deriv, double compensateBias);
    /// MuJoCo model
    const mjModel *m;
    /// Perturbation for central differences
    double eps = 1e-6;
    /// Variables for differentiation
    eigVd xNewTemp, xNewP, xNewN;
    eigVd uTemp;
    /// Objects
    Params *cp;
    Control *cc;
};

#endif //NUMDIFF_H
