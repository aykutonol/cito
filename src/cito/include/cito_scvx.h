// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_SCVX class consists of functions that setup the successive
// convexification algorithm.

#ifndef CITO_SCVX_H
#define CITO_SCVX_H

class CitoNumDiff
{
private:
    int maxIter = 25;
    double Jtemp[maxIter+1], J[maxIter+1], L[maxIter+1];
    double r[maxIter+1], rho[maxIter+1], dL[maxIter+1], dJ[maxIter+1];
    r[0] = 1e2;   // initial trust region radius
    double dLTol = 1e-4;
    double rho0 = 0, rho1 = 0.25, rho2 = 0.90, rMin = 0, rMax = 1e20;
    double alpha = 2, beta = 3.2;
    int accept[maxIter+1], dLTolMet = 0, stop = 0;

public:
    CitoScvx();
    ~CitoScvx();
};

#endif //CITO_SCVX_H
