/*! Parameters */
/**
 *  \brief Params contains definitions that are used across CITO classes.
 *
 *  This header contains global types, structs, and paths as well as the
 *  Params class that parses the model and config files and defines
 *  parameters and utility functions that are used across classes.
 *
 *  \author Aykut Onol
 */

#ifndef PARAMS_H
#define PARAMS_H

#include <iostream>
#include <string>
#include <chrono>
#include <ctime>

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

/// Params class
class Params
{
public:
    /// Constructor
    Params(const mjModel *model_);
    /// Destructor
    ~Params();
    /// Simulation and model parameters
    const mjModel *model;
    double tf, tc, dt;
    int N, ndpc, nu, nv, n, m, nTraj, nPair,
        nFree, *pFree, *bFree, *dAct, *pAct,
        *quatAdr, *dofAdr,
        taskType;
    std::vector<Eigen::Vector2i> sites;
    eigMd nCS0;
    mjtNum unit_x[3] = {1., 0., 0.}; // unit-x vector
    /// Task parameters
    eigVd desiredPos;
    int controlJointDOF0;
    double weight[4];
    /// Utility functions
    /// This function returns the skew symmetric matrix representation of a 3D vector
    Eigen::Matrix3d skew(const Eigen::Vector3d &a);
    /// This function performs and returns a x b using skew-symmetric transformation
    Eigen::Vector3d skewCross(const Eigen::Vector3d &a, const Eigen::Vector3d &b);
    /// This function converts a quaternion (w,x,y,z) into Euler angles
    Eigen::Vector3d quat2Euler(const Eigen::Vector4d &q);
    /// This function calculates the contact normal Jacobian w.r.t. rotational DOF
    Eigen::Matrix3d evalNormalJac(const Eigen::Vector4d &q, int pair);
    /** This function converts free joints' quaternions to Euler angles so that
     *  the dimensionality of the state vector is 2*nv instead of nq+nv */
    eigVd getState(const mjData *d);
};

#endif //PARAMS_H
