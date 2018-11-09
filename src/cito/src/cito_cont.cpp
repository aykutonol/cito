// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONT class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

// ***** CLASS TYPE ************************************************************
// Robot and environment specific

#include "cito_cont.h"

// intermediate control variables **********************************************


//setControl: sets generalized forces on joints and bodies
void SCvx::setControl(const mjModel* m, mjData* d, const ctrlVec_t u)
{
    // set control
    for( int i=0; i<ndof; i++ )
    {
      d->ctrl[i] = u[i] + d->qfrc_bias[u_id[i]];
    }
    // get kcon
    double kcon = u[ndof];
    // contact model
    ho = this->contactModel(m, d, kcon);
    // set external force on the torso
    for( int i=0; i<6; i++ )
    {
      d->xfrc_applied[mj_name2id(m, mjOBJ_BODY, "torso")*6+i] = ho[i];
    }
}

// contactModel: returns contact wrench given kcon
Eigen::Matrix<double, 6, 1> SCvx::contactModel(const mjModel* m, const mjData* d, double kcon)
{
    for( int c=0; c<ncc; c++ )
    {
        // contact normal
        for( int i=0; i<4; i++ )
        {
          obj_q[i] = d->qpos[i+3];
        }
        mju_rotVecQuat(cc_n, cc_n0, obj_q);
        // distance
        for( int i=0; i<3; i++ )
        {
          pcc(i) = d->site_xpos[mj_name2id(m, mjOBJ_SITE, "site_cc1")*3+i];
          pee(i) = d->site_xpos[mj_name2id(m, mjOBJ_SITE, "site_ee")*3+i];
          cc_n_Eig(i) = cc_n[i];
        }
        pee2cc = pcc-pee;
        phi_e = sqrt(pow(pee2cc[0],2)+pow(pee2cc[1],2)+pow(pee2cc[2],2));
        phi_n = pee2cc.dot(cc_n_Eig);
        zeta = tanh(200*phi_e);
        phi_c = zeta*phi_e + (1-zeta)*phi_n;
        // normal contact force
        fn = kcon*exp(-acon*phi_c);
        // object wrench & joint torques
        lambda = fn*cc_n_Eig;
        // wrench on the object
        for( int i=0; i<3; i++ )
        po(i) = d->qpos[i];
        L = pee - po;
        Lhat <<     0,   -L[2],    L[1],
                 L[2],       0,   -L[0],
                -L[1],    L[0],       0;
        W << Eigen::MatrixXd::Identity(3,3),
             Lhat;
    }

    return wrench;
}
