/*! Parameters */
/**
 *  \brief CitoParams contains definitions that are used across classes.
 *
 *  This header contains global types, structs, and paths as well as the
 *  CitoParams class that parses the model and config files and defines
 *  parameters that are used across classes.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_PARAMS_H
#define CITO_PARAMS_H

#include <iostream>
#include <chrono>
#include <string>

#include "mujoco.h"
#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>

/// Types
typedef Eigen::VectorXd eigVd;
typedef Eigen::MatrixXd eigMd;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, 1> eigVm;
typedef Eigen::Matrix<mjtNum, Eigen::Dynamic, Eigen::Dynamic> eigMm;
typedef std::vector<eigMm> eigTm;

/// Structs
struct trajectory
{
    eigMm X;
    eigMd U;
    eigTm Fx;
    eigTm Fu;
};

/// Paths
namespace paths
{
    const std::string workspaceDir = std::getenv("CITO_WS");
}

/// CitoParams class
class CitoParams{
public:
    /// Constructor
    CitoParams(const mjModel *model);
    /// Destructor
    ~CitoParams() {}
    /// Simulation and model parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct,
        *quatAdr, *dofAdr;
    Eigen::VectorXi sPair1, sPair2;
    eigMm nCS;
    /// Task parameters
    eigVd desiredPos;
    int controlJointDOF0;
    double weight[4];
};

#endif //CITO_PARAMS_H
