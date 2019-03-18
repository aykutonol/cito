/*! Control and Contact Model */
/**
 *  \brief CitoControl class consists of functions for control and contact model
 *
 *  This class defines functions for calculating and setting control- and contact-
 *  related variables, i.e., joint torques and external forces on the free bodies.
 *
 *  \author Aykut Onol
 */

#include "cito_params.h"
#include "mj_savelog.h"

#ifndef CITO_CONTROL_H
#define CITO_CONTROL_H

class CitoControl
{
public:
    /// Constructor
    CitoControl(const mjModel* model);
    /// Destructor
    ~CitoControl();
    /// This function takes a full control step given a control input
    void takeStep(mjData*d, const eigDbl u, bool save);
    /// This function sets generalized forces on joints and free bodies
    void setControl(mjData* d, const eigDbl u);
    /** This function converts free joints' quaternions to Euler angles so that
     *  the dimensionality of the state vector is 2*nv instead of nq+nv */
    eigMjc getState(const mjData* d);
    /// This function gets bounds on joint positions, actuator forces from the model
    void getBounds();
    /// position & torque limits
    double *qposLB, *qposUB, *tauLB, *tauUB;
    int    *isJFree, *isAFree;

private:
    /// This function returns contact wrench given current state and control input
    eigDbl contactModel(const mjData* d, const eigDbl u);
    /// Contact wrench
    eigDbl h, hCon;
    /// Contact model variables
    double phiE, phiN, zeta, phiC, gamma, alpha, phiR;
    Eigen::Matrix<double, 3, 1> pSR, pSE, pBF, nCS, vRE, vEF, lambda;
    /// getState variables
    eigMjc x;
    Eigen::Matrix<mjtNum, 4, 1> jFreeQuat;
    /// Objects to be initialized
    const mjModel* m;
    MjSaveLog  sl;
    CitoParams cp;
};

#endif //CITO_CONTROL_H
