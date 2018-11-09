// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //
#include "cito_params.h"

#ifndef CITO_TYPES_H
#define CITO_TYPES_H

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

#endif //CITO_TYPES_H
