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
private:
    // SCvx parameters
    int maxIter = 25;  // maximum number of iterations
    double r0 = 1e2;   // initial trust region radius
    double Jtemp[maxIter+1], J[maxIter+1], L[maxIter+1];
    double r[maxIter+1], rho[maxIter+1], dL[maxIter+1], dJ[maxIter+1];
    double dLTol = 1e-4;
    double rho0 = 0, rho1 = 0.25, rho2 = 0.90, rMin = 0, rMax = 1e20;
    double alpha = 2, beta = 3.2;
    bool accept[maxIter+1], dLTolMet = 0, stop = 0;
    // state/control vectors/matrices
    stateVecThread X, dX, XL;   ctrlVecThread  U, dU, Utemp;
    stateMatThread Fx;          ctrlMatThread  Fu;
    // getCost
    Eigen::VectorXd po_f(6);
    Eigen::Matrix<double, NTS, params::npair> k;
    double Jt, Ji, J;
    // custom objects
    CitoControl cc;
    CitoNumDiff nd;

public:
    CitoSCvx();
    ~CitoSCvx() {}
};

#endif //CITO_SCVX_H
