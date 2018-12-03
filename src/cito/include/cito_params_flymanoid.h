// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_PARAMS_XYZ class defines parameters and types that are customized with
// respect to the robot and the environment.

// ***** CLASS TYPE ************************************************************
// Robot and environment specific

// SPECIFIC to flymanoid.xml

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "mujoco.h"
#include <Eigen/Dense>

#ifndef CITO_PARAMS_FLYMANOID_H
#define CITO_PARAMS_FLYMANOID_H


// PARAMS **********************************************************************
const double tf = 2.00;           // [s] final time
const double tc = 1e-1;           // [s] control sampling period
const double dt = 2e-3;           // [s] dynamic sampling period
const int ncts  = floor(tf/tc);   // number of control time steps
const int ndts  = floor(tf/dt)+1; // number of dynamic time steps
const int ndpc  = floor(tc/dt);   // number of dynamic time steps per a control step
const int nact  = 8;              // number of actuated dof
const int nfree = 1;              // number of free joints

// contact related
const int ncrbt = 4;              // number of contact candidates on the robot (end effectors)
const int ncenv = 4;              // number of contact candidates in the environment
const int nsite = ncrbt + ncenv;  // number of sites
const int npair = 8;              // number of contact pairs (4 for each side)

// specific indices
const int jact[nact]    = {6,7,8,9,10,11,12,13};  // actuated dof indices
const int jfree[nfree]  = {0};                    // index of the free body
const int spair1[npair] = {4,4,5,5,6,6,7,7};      // indices of the sites on the robot
const int spair2[npair] = {1,3,0,2,1,3,0,2};      // corresponding indices of the sites in the environment

// initial pose of the robot and controls
Eigen::Matrix<mjtNum, nact, 1> qpos0;
double kcon0[npair] = {10,10,10,10,10,10,10,10};
double acon[npair]  = {20,20,20,20,20,20,20,20};

// contact surface normals for each pair
double csn[npair*3] = {0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0};

// constant variables for later use
const int NV  = nact + 6*nfree;
const int N   = 2*NV;
const int M   = nact + npair;
const int NTS = ncts;

// TYPES ***********************************************************************
// instantaneous eigen+mujoco types
typedef Eigen::Matrix<mjtNum, NV, NV, Eigen::RowMajor>  dofMat_t;;
typedef Eigen::Matrix<mjtNum, N, N,   Eigen::RowMajor>  stateMat_t;
typedef Eigen::Matrix<mjtNum, N, 1>                     stateVec_t;
typedef Eigen::Matrix<mjtNum, NV, M,  Eigen::RowMajor>  dofCtrlMat_t;
typedef Eigen::Matrix<mjtNum, N, M,   Eigen::RowMajor>  ctrlMat_t;
typedef Eigen::Matrix<mjtNum, M, 1>                     ctrlVec_t;


// threaded types for multiple time steps
typedef std::vector<stateVec_t, Eigen::aligned_allocator<stateVec_t>> stateVecThread;
typedef std::vector<stateMat_t, Eigen::aligned_allocator<stateMat_t>> stateMatThread;
typedef std::vector<ctrlVec_t,  Eigen::aligned_allocator<ctrlVec_t>>  ctrlVecThread;
typedef std::vector<ctrlMat_t,  Eigen::aligned_allocator<ctrlMat_t>>  ctrlMatThread;


#endif //CITO_PARAMS_FLYMANOID_H
