/*! Control and Contact Model */
/**
 *  \brief Control consists of methods for control and contact model
 *
 *  This class defines functions for calculating and setting control- and contact-
 *  related variables, e.g., joint torques and external forces on the free bodies.
 *
 *  \author Aykut Onol
 */

#ifndef CONTROL_H
#define CONTROL_H

#include "cito/params.h"
#include "cito/savelog.h"

#include <fcl/fcl.h>

class Control
{
public:
    /// Constructor
    Control(const mjModel *m_, Params *cp_);
    /// Destructor
    ~Control();
    /// This function returns the 3D transform of a desired site
    fcl::Transform3d getSiteTransform(const mjData *d, int site_id);
    /// This function creates a collision geometry given a site ID
    std::shared_ptr<fcl::CollisionGeometryd> createCollGeom(const mjModel *m, int site_id);
    /// This function returns FCL distance calculation results for all contact pairs
    std::vector<fcl::DistanceResultd> calcDistance(const mjData *d);
    /// This function takes a full control step given a control input
    void takeStep(mjData *d, const eigVd &u, int save, double compensateBias);
    /// This function sets generalized forces on joints and free bodies
    void setControl(mjData *d, const eigVd &u, double compensateBias);
    /// This function gets bounds on joint positions, actuator forces from the model
    void getBounds();
    /// Position & torque limits
    double *qposLB, *qposUB, *tauLB, *tauUB;
    int *isJFree, *isAFree;
    // Contact model curvature
    double alpha;

private:
    /// This function returns contact wrench given current state and control input
    eigVd contactModel(const mjData *d, const eigVd &u);
    /// MuJoCo model
    const mjModel *m;
    /// Contact wrench
    eigVd h, hCon;
    /// Contact model variables
    double gamma;
    Eigen::Vector3d nCS, lambda, pCoM, r, sEnvPos;
    /// Objects
    Params *cp;
    SaveLog sl;
    /// FCL solver and variables
    fcl::DistanceRequestd distReq;
    fcl::DistanceResultd distRes;
    fcl::detail::GJKSolver_libccdd fclDist;
    std::unordered_map<int, fcl::CollisionObjectd *> collObjs;
};

#endif //CONTROL_H
