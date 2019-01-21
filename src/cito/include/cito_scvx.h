// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SCVX class consists of functions that setup the successive
// convexification algorithm.

#ifndef CITO_SCVX_H
#define CITO_SCVX_H

#include <iostream>
#include "cito_numdiff.h"
#include "cito_sqopt.h"

class CitoSCvx
{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
    CitoSCvx(const mjModel* model);
    ~CitoSCvx() {}
    // ***** FUNCTIONS *************************************************************
    double getCost(const stateVec_t Xfinal, const ctrlVecThread U);
    void runSimulation(const ctrlVecThread U0, bool linearize, bool save);
    void solveSCvx(const ctrlVecThread U);

private:
    // ***** PARAMETERS ************************************************************
    const mjModel* m;
    // SCvx parameters
    static const int maxIter = 2;  // maximum number of iterations
    double r0 = 1e2;                // initial trust region radius
    double Jtemp[maxIter+1], J[maxIter+1], L[maxIter+1];
    double r[maxIter+1], rho[maxIter+1], dL[maxIter+1], dJ[maxIter+1];
    double dLTol = 1e-4;
    double rho0 = 0, rho1 = 0.25, rho2 = 0.90, rMin = 0, rMax = 1e20;
    double alpha = 2, beta = 3.2;
    bool accept[maxIter+1], dLTolMet = 0, stop = 0;
    // state/control vectors/matrices for the succession, change, and linear approximation
    stateVecThread Xs, dX, Xl;      ctrlVecThread  Us, dU, Utemp;
    stateMatThread Fx;              ctrlMatThread  Fu;
    // getCost
    Eigen::Matrix<double, 6, 1> finalPose;
    kconVecThread Kcon;
    double KconSN;
    double Jf, Ji, Jt;      // final, integrated, and total cost values
    // ***** OBJECTS ***************************************************************
    CitoControl cc;
    CitoNumDiff nd;
    CitoSQOPT   sq;
};

#endif //CITO_SCVX_H
