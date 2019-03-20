/*! Parameters */
/**
 *  \brief CitoParams contains parameter definitions
 *
 *  This class parses the model and config files and defines parameters
 *  as well as types and structs that are used across classes.
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

/// Types
typedef Eigen::VectorXd eigVd;
typedef Eigen::MatrixXd eigMd;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, Eigen::Dynamic> eigMm;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> derTraj;

/// Structs
struct trajectory
{
    eigMm X;
    eigMd  U;
    derTraj Fx;
    derTraj Fu;
};

/// CitoParams class
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
    eigMm nCS;
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

#endif //CITO_PARAMS_H
