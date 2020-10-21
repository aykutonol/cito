/*! Successive Convexification */
/**
 *  \brief SCVX consists of methods to run the SCVX algorithm.
 *
 *  This class defines functions that are used to roll-out the dynamics,
 *  evaluate the cost, and execute the SCVX algorithm.
 *
 *  \author Aykut Onol
 */

#ifndef SCVX_H
#define SCVX_H

#include "cito/numdiff.h"
#include "cito/sqopt.h"

class SCVX
{
public:
    /// Constructor
    SCVX(const mjModel *m_, Params *cp_, Control *cc_);
    /// Destructor
    ~SCVX();
    /// This function returns the nonlinear cost given control trajectory and final state
    double getCost(const eigMd &X, const eigMd &U);
    /// This function rolls-out and linearizes the dynamics given control trajectory
    trajectory runSimulation(const eigMd &U0, bool linearize, bool save, double compensateBias);
    /// This function executes the successive convexification algorithm
    eigMd solveSCVX(const eigMd &U);
    // This function refreshes SCVX variables for a new run
    void refresh();

private:
    /// MuJoCo model
    const mjModel *m;
    /// SCVX parameters
    int maxIter;                  // maximum number of iterations
    double *J, *JTemp, *JTilde,   // cost terms
        *r, *dJ, *dL, *rho,       // trust region radius, change, and similarity
        dLTol,                    // stopping criteria in terms of dL
        rho0, rho1, rho2,         // similarity thresholds
        beta_expand, beta_shrink, // trust-region expand and shrink factors
        r0, rMin, rMax;           // initial, min, and max trust-region radius
    bool *accept, dLTolMet = false, stop = false;
    /// Trajectories
    eigMd XSucc, dX, XTilde;
    eigMd USucc, dU, UTemp;
    eigTd Fx, Fu;
    trajectory traj, trajS, trajTemp;
    /// Cost function variables
    eigVd finalPos;
    double Jf, Ji, Jt; // final, integrated, and total cost values
    /// Objects
    Params *cp;
    Control *cc;
    NumDiff nd;
    SQOPT sq;
};

#endif //SCVX_H
