// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_PARAMS class defines parameters and types that are customized w.r.t.
// the robot and the environment.

// ***** CLASS TYPE ************************************************************
// Robot and environment specific

// This file is SPECIFIC to the model flymanoid.xml

#include "mujoco.h"
#include <Eigen/Dense>

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

// properties of the joint(s) to be controlled
const double desiredPose[6] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const int controlJointPos0 = 0;   // index of the first element of the joint position

// PARAMS **********************************************************************
namespace params
{
  const double tf = 0.20;           // [s] final time
  const double tc = 1e-1;           // [s] control sampling period
  const double dt = 2e-3;           // [s] dynamic sampling period
  const int ncts  = floor(tf/tc);   // number of control time steps
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
  const int jfree[nfree]  = {0};                    // index of the free joint
  const int bfree[nfree]  = {5};                    // index of the free body
  const int spair1[npair] = {4,4,5,5,6,6,7,7};      // indices of the sites on the robot
  const int spair2[npair] = {1,3,0,2,1,3,0,2};      // corresponding indices of the sites in the environment

  // initial pose of the robot and controls
  const double kcon0[npair] = {10,10,10,10,10,10,10,10};
  const double acon[npair]  = {10,10,10,10,10,10,10,10};
  const double phi_r = 200;     // parameter for the radius of the distance sphere (200 ~ 1 cm)

  // contact surface normals for each pair
  const double csn[npair*3] = {0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0, 0,1,0};
}

// constant variables for types
const int NU    = params::nact;           // number of actuated joints
const int NPAIR = params::npair;          // number of contact pairs
const int NV    = NU + 6*params::nfree;   // degrees of freedom
const int N     = 2*NV;                   // dimensionality of states
const int M     = NU + NPAIR;             // dimensionality of controls
const int NTS   = params::ncts;           // number of control time
const int NTRAJ = (NTS+1)*N + NTS*M;      // number of trajectory variables

// TYPES ***********************************************************************
// instantaneous eigen+mujoco types
typedef Eigen::Matrix<mjtNum, N, N, Eigen::ColMajor> stateMat_t;
typedef Eigen::Matrix<mjtNum, N, 1>                  stateVec_t;
typedef Eigen::Matrix<mjtNum, N, M, Eigen::ColMajor> ctrlMat_t;
typedef Eigen::Matrix<mjtNum, M, 1>                  ctrlVec_t;
typedef Eigen::Matrix<mjtNum, NPAIR, 1>              kconVec_t;
// threaded types for multiple time steps
typedef std::vector<stateVec_t, Eigen::aligned_allocator<stateVec_t>> stateVecThread;
typedef std::vector<stateMat_t, Eigen::aligned_allocator<stateMat_t>> stateMatThread;
typedef std::vector<ctrlVec_t,  Eigen::aligned_allocator<ctrlVec_t>>  ctrlVecThread;
typedef std::vector<ctrlMat_t,  Eigen::aligned_allocator<ctrlMat_t>>  ctrlMatThread;
typedef std::vector<kconVec_t,  Eigen::aligned_allocator<kconVec_t>>  kconVecThread;


#endif //CITO_PARAMS_H
