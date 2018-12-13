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
private:
    const mjModel* m;
    // control variables
    Eigen::Matrix<double, params::npair, 1>   kcon;
    Eigen::Matrix<double, 6*params::nfree, 1> hcon;
    // contact model variables
    double phi_e, phi_n, zeta, phi_c, fn;
    Eigen::Matrix<double, 3, 1> p_sr, p_se, p_bf, n_cs, v_re, v_ef, lambda;
    Eigen::Matrix<double, 6*params::nfree, 1> h;
    // variables for getState
    stateVec_t x;
    Eigen::Matrix<mjtNum, 4, 1> obj_q;
    // internal functions
    Eigen::Matrix<double, 6*params::nfree, 1> contactModel(const mjData* d, const ctrlVec_t u);

public:
    CitoControl(const mjModel* model);
    ~CitoControl() {}
    void setControl(mjData* d, const ctrlVec_t u);
    stateVec_t getState(const mjData* d);
};

#endif //CITO_CONTROL_H
