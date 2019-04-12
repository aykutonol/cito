/*! Successive Convexification */
/**
 *  \brief CitoSCvx class consists of functions to run the SCvx algorithm
 *
 *  This class defines functions that are used to roll-out the dynamics,
 *  evaluate the cost, and execute the SCvx algorithm.
 *
 *  \author Aykut Onol
 */


#ifndef CITO_SCVX_H
#define CITO_SCVX_H

#include "cito_sqopt.h"
#include "cito_control.h"
#include "cito_deriv.h"


class CitoSCvx
{
public:
    /// Constructor
    CitoSCvx(const mjModel* model);
    /// Destructor
    ~CitoSCvx() {}
    /// This function returns the nonlinear cost given control trajectory and final state
    double getCost(const eigMm X, const eigMd U);
    /// This function rolls-out and linearizes the dynamics given control trajectory
    trajectory runSimulation(const eigMd U0, bool linearize, bool save, double compensateBias);
    /// This function executes the successive convexification algorithm
    eigMd solveSCvx(const eigMd U);

private:
    /// MuJoCo model
    const mjModel* m;
    /// SCvx parameters
    int maxIter;                        // maximum number of iterations
    double *J, *JTemp, *JTilde,         // cost terms
           *r, *dJ, *dL, *rho,          // trust region radius, change, and similarity
           dLTol,                       // stopping criteria in terms of dL
           rho0, rho1, rho2,            // similarity thresholds
           beta_expand, beta_shrink,    // trust-region expand and shrink factors
           rMin, rMax;                  // trust-region radius limits
    bool *accept, dLTolMet = false, stop = false;
    /// Trajectories
    eigMm XSucc, dX, XTilde;
    eigMd USucc, dU, UTemp;
    eigTm Fx, Fu;
    trajectory traj, trajS, trajTemp;
    /// Cost function variables
    eigVd finalPos;
    double Jf, Ji, Jt;      // final, integrated, and total cost values
    /// Objects
    CitoParams  cp;
    CitoControl cc;
    CitoDeriv   cd;
    CitoSQOPT   sq;
};

#endif //CITO_SCVX_H
