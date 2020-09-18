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
    eigMd nCS0;
    /// Task parameters
    eigVd desiredPos;
    int controlJointDOF0;
    double weight[4];
    /// Utility functions
    /// This function returns the skew symmetric matrix representation of a 3D vector
    Eigen::Matrix3d skew(const Eigen::Vector3d& a);
    /// This function performs and returns a x b using skew-symmetric transformation
    Eigen::Vector3d skewCross(const Eigen::Vector3d& a, const Eigen::Vector3d& b);
    /// This function converts a quaternion (w,x,y,z) into Euler angles
    Eigen::Vector3d quat2Euler(const Eigen::Vector4d& q);
    /// This function calculates the contact normal Jacobian w.r.t. rotational DOF
    Eigen::Matrix3d evalNormalJac(const Eigen::Vector4d& q, int pair);
};

#endif //CITO_PARAMS_H
