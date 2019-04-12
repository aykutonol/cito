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
    CitoDeriv(const mjModel* model);
    /// Destructor
    ~CitoDeriv() {}
    /// This function calculates the derivatives of the dynamics given state and control vectors
    void linDyn(const mjData* dMain, const eigVd uMain, mjtNum* Fxd, mjtNum* Fud, double compensateBias);
private:
    /// MuJoCo model
    const mjModel* m;
    /// Objects
    CitoParams  cp;
};


#endif //CITO_DERIV_H
