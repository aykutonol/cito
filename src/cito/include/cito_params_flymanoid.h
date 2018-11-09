// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// SPECIFIC to flymanoid.xml

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "mujoco.h"
#include <Eigen/Dense>

#ifndef CITO_PARAMS_FLYMANOID_H
#define CITO_PARAMS_FLYMANOID_H

// initialization
const double tf = 2.00;           // [s] final time
const double tc = 1e-1;           // [s] control sampling period
const double dt = 2e-3;           // [s] dynamic sampling period
const int ncts  = floor(tf/tc);   // number of control time steps
const int ndts  = floor(tf/dt)+1; // number of dynamic time steps
const int ndpc  = floor(tc/dt);   // number of dynamic time steps per a control step
const int ndof  = 8;              // number of actuated dof
const int nfree = 1;              // number of free joints
const int ncc   = 8;              // number of contact pairs

// controlled dof indices
const int u_id[ndof] = {6,7,8,9,10,11,12,13};

// initial pose of the robot and controls
Eigen::Matrix<mjtNum, ndof, 1> qpos0;
double kcon0 = 10;
double acon  = 20;

// sawyer tabletop (1 object) model
const int NV  = ndof + 6*nfree;
const int N   = 2*NV;
const int M   = ndof+ncc;
const int NTS = ncts;

#endif //CITO_PARAMS_FLYMANOID_H
