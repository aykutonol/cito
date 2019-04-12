/*! Analytical Derivatives */
/**
 *  \brief CitoDeriv class consists of functions for taking derivatives.
 *
 *  This class defines functions for calculating the derivatives of the
 *  dynamics using Pinocchio.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_DERIV_H
#define CITO_DERIV_H

#include "cito_params.h"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"


class CitoDeriv {
public:
    /// Constructor
    CitoDeriv(const mjModel* mModel);
    /// Destructor
    ~CitoDeriv() {}
    /// This function calculates the derivatives of the dynamics given state and control vectors
    void linDyn(const mjData* d, const eigVd u, mjtNum* Fxd, mjtNum* Fud, double compensateBias);
private:
    /// MuJoCo model
    const mjModel* m;
    /// Pinocchio model, state, and control
    pinocchio::Model model;
    eigVd q, v, tau;
    /// Derivative matrices
    eigMm Fx, Fu;
    /// Time coefficient: (time step size)*(number of dynamic time steps per control period)
    double tM;
    /// Objects
    CitoParams cp;
};

#endif //CITO_DERIV_H
