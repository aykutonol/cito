// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONTROL class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

#include "cito_params.h"

#ifndef CITO_CONTROL_H
#define CITO_CONTROL_H

class CitoControl
{
public:
    // ***** CONSTRUCTOR/DESTRUCTOR ************************************************
    CitoControl(const mjModel* model);
    ~CitoControl();
    // ***** FUNCTIONS *************************************************************
    void takeStep(mjData*d, const ctrlVec_t u);
    void setControl(mjData* d, const ctrlVec_t u);
    stateVec_t getState(const mjData* d);
    void getBounds();
    // ***** PARAMETERS ************************************************************
    // position & torque limits
    double *qposLB  = new double[NV];
    double *qposUB  = new double[NV];
    double *tauLB   = new double[NU];
    double *tauUB   = new double[NU];
    int    *isJFree = new int[NV];
    int    *isAFree = new int[NU];

private:
    // ***** FUNCTIONS *************************************************************
    Eigen::Matrix<double, 6*params::nfree, 1> contactModel(const mjData* d, const ctrlVec_t u);
    // ***** PARAMETERS ************************************************************
    const mjModel* m;
    // control variables
    Eigen::Matrix<double, 6*params::nfree, 1> h, hCon;
    // contact model variables
    double phiE, phiN, zeta, phiC, gamma;
    Eigen::Matrix<double, 3, 1> pSR, pSE, pBF, nCS, vRE, vEF, lambda;
    // variables for getState
    stateVec_t x;
    Eigen::Matrix<mjtNum, 4, 1> jfree_quat;
};

#endif //CITO_CONTROL_H
