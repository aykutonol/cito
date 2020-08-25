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
#include <string>
#include <chrono>

#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>
#include "mujoco.h"

/// Types
typedef Eigen::VectorXd eigVd;
typedef Eigen::MatrixXd eigMd;
typedef std::vector<eigMd> eigTd;

/// Structs
struct trajectory
{
    eigMd X;
    eigMd U;
    eigTd Fx;
    eigTd Fu;
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
    CitoParams(const mjModel *model_);
    /// Destructor
    ~CitoParams();
    /// Simulation and model parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct,
        *quatAdr, *dofAdr;
    std::vector<Eigen::Vector2i> sites;
    eigMd nCS;
    /// Task parameters
    eigVd desiredPos;
    int controlJointDOF0;
    double weight[4];
    /// Utility functions
    Eigen::Vector3d skewCross(const Eigen::Vector3d& a, const Eigen::Vector3d& b);
};

#endif //CITO_PARAMS_H
