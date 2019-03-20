/*! Parameters */
/**
 *  \brief CitoParams contains user-specific and general definitions
 *
 *  This class defines global variables that are specific to simulation,
 *  robot, and environment as well as general types and structures.
 *
 *  \author Aykut Onol
 */

#include <stdio.h>
#include <iostream>
#include <string>

#include "mujoco.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <yaml-cpp/yaml.h>

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

typedef Eigen::MatrixXd eigDbl;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, Eigen::Dynamic> eigMjc;

class CitoParams{
public:
    /// Constructor
    CitoParams(const mjModel *model);
    /// Destructor
    ~CitoParams() {}
    /// Parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct,
        *quatAdr, *dofAdr;
    Eigen::VectorXi sPair1, sPair2;
    eigMjc nCS;
};

/// User-specific paths
namespace paths {
    const std::string workspaceDir = std::getenv("CITO_WS");
    const std::string modelFile    = "sawyer_push.xml";
}

/// Specify the parameters related to simulation, model, and contact candidates
/// Parameters
namespace params {
    /// Model
    const int nact  = 7;                    // number of actuated joints
    const int nfree = 1;                    // number of free joints
    /// Contact
    const int npair = 1;                // number of contact pairs
}
// The following constants and types are not changed */
/// Constants
const int NU    = params::nact;         // number of actuated joints
const int NPAIR = params::npair;        // number of contact pairs
const int NV    = NU + 6*params::nfree; // degrees of freedom
const int N     = 2*NV;                 // dimensionality of states
const int M     = NU + NPAIR;           // dimensionality of controls
/// Types
// eigen+mujoco types for a time instant
typedef Eigen::Matrix<mjtNum, N, N, Eigen::ColMajor> stateDer;
typedef Eigen::Matrix<mjtNum, N, M, Eigen::ColMajor> ctrlDer;
// threaded types for trajectories
typedef std::vector<stateDer, Eigen::aligned_allocator<stateDer>> stateDerTraj;
typedef std::vector<ctrlDer,  Eigen::aligned_allocator<ctrlDer>>  ctrlDerTraj;
/// Structs
struct trajectory
{
    eigMjc X;
    eigDbl  U;
    stateDerTraj Fx;
    ctrlDerTraj  Fu;
};

#endif //CITO_PARAMS_H
