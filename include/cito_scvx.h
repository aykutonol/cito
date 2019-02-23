// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SCVX class consists of functions that setup the successive
// convexification algorithm.

#ifndef CITO_SCVX_H
#define CITO_SCVX_H

#include <iostream>
#include <chrono>
#include "cito_numdiff.h"
#include "cito_sqopt.h"


class CitoSCvx
{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
    CitoSCvx(const mjModel* model);
    ~CitoSCvx() {}
    // ***** FUNCTIONS *************************************************************
    double getCost(const stateVec_t XFinal, const ctrlVecThread U);
    trajectory runSimulation(const ctrlVecThread U0, bool linearize, bool save);
    ctrlVecThread solveSCvx(const ctrlVecThread U);

private:
    // ***** PARAMETERS ************************************************************
    const mjModel* m;
    // SCvx parameters
    int maxIter;                        // maximum number of iterations
    double *J, *JTemp, *JTilde,         // cost terms
           *r, *dJ, *dL, *rho,          // trust region radius, change, and similarity
           dLTol,                       // stopping criteria in terms of dL
           rho0, rho1, rho2,            // similarity thresholds
           beta_expand, beta_shrink,    // trust-region expand and shrink factors
           rMin, rMax;                  // trust-region radius limits
    bool *accept, dLTolMet = false, stop = false;
    // state/control vectors/matrices for the succession, change, and linear approximation
    stateVecThread XSucc, dX, XTilde;   ctrlVecThread USucc, dU, UTemp;
    stateDerThread Fx;                  ctrlDerThread Fu;
    // trajectories
    trajectory traj, trajS, trajTemp;
    // getCost
    double weight[4];
    int controlJointDOF0;
    Eigen::Matrix<double, 6, 1> desiredPose, desiredVelo, finalPose, finalVelo;
    kConVecThread KCon;
    double KConSN;          // total squared norm of virtual stiffness variables
    double Jf, Ji, Jt;      // final, integrated, and total cost values
    // ***** OBJECTS ***************************************************************
    CitoControl cc;
    CitoNumDiff nd;
    CitoSQOPT   sq;
};

#endif //CITO_SCVX_H
