/*! Control and Contact Model */
/**
 *  \brief CitoControl class consists of functions for control and contact model
 *
 *  This class defines functions for calculating and setting control- and contact-
 *  related variables, i.e., joint torques and external forces on the free bodies.
 *
 *  \author Aykut Onol
 */

#ifndef CITO_CONTROL_H
#define CITO_CONTROL_H

#include "cito_params.h"
#include "mj_savelog.h"

#include <fcl/fcl.h>

class CitoControl
{
public:
    /// Constructor
    CitoControl(const mjModel* m_, CitoParams* cp_);
    /// Destructor
    ~CitoControl();
    /// This function takes a full control step given a control input
    void takeStep(mjData* d, const eigVd u, bool save, double compensateBias);
    /// This function sets generalized forces on joints and free bodies
    void setControl(mjData* d, const eigVd u, double compensateBias);
    /** This function converts free joints' quaternions to Euler angles so that
     *  the dimensionality of the state vector is 2*nv instead of nq+nv */
    eigVd getState(const mjData* d);
    /// This function gets bounds on joint positions, actuator forces from the model
    void getBounds();
    /// Position & torque limits
    double *qposLB, *qposUB, *tauLB, *tauUB;
    int    *isJFree, *isAFree;
    /// FCL functions
    fcl::Transform3d getSiteTransform(const mjData* d, int site_id);
    std::shared_ptr<fcl::CollisionGeometryd> createCollGeom(const mjModel* m, int site_id);

private:
    /// This function returns contact wrench given current state and control input
    eigMd contactModel(const mjData* d, const eigVd u);
    /// MuJoCo model
    const mjModel* m;
    /// Contact wrench
    eigMd h, hCon;
    /// Contact model variables
    double phiE, phiN, zeta, phiC, gamma, alpha, phiR;
    Eigen::Vector3d pSR, pSE, pBF, nCS, vRE, vEF, lambda;
    mjtNum unit_x[3] = {1., 0., 0.};    // unit-x vector
    /// getState variables
    eigVd x;
    Eigen::Matrix<mjtNum, 4, 1> jFreeQuat;
    /// Objects
    CitoParams *cp;
    MjSaveLog sl;
    /// FCL solver and variables
    fcl::DistanceRequestd distReq;
    fcl::DistanceResultd  distRes;
    fcl::detail::GJKSolver_libccdd fclDist;
    std::unordered_map<int, fcl::CollisionObjectd*> collObjs;
};

#endif //CITO_CONTROL_H
