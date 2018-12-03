// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONT class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

// ***** CLASS TYPE ************************************************************
// Robot, environment, and physics engine specific

#include "cito_cont.h"

//setControl: sets generalized forces on joints and free bodies
void CitoCont::setControl(mjData* d, const ctrlVec_t u)
{
    // set control given the control input
    for( int i=0; i<nact; i++ )
    {
      d->ctrl[i] = u[i] + d->qfrc_bias[act_id[i]];
    }
    // get kcon from the control input
    for( int i=0; i<npair; i++ )
    {
      kcon[i] = u[nact+i];
    }
    // contact model
    hcon[6*nfree] = this->contactModel(d, kcon);
    // set external forces on the free bodies
    for( int i=0; i<nfree; i++ )
    {
      for( int j=0; j<6; j++ )
      {
        d->xfrc_applied[free_id[i]*6+j] = hcon[i*6+j];
      }
    }
}

// contactModel: returns contact wrench given kcon
Eigen::Matrix<double, 6*nfree, 1> CitoCont::contactModel(const mjData* d, double kcon)
{
    h.setZero();
    // loop for each contact pair
    for( int p_i=0; p_i<npair; p_i++ )
    {
        // vectors in the world frame
        for( int i=0; i<3; i++ )
        {
          p_sr[i] = d->site_xpos[spair1[p_i]*3+i];  // position of the site on the robot
          p_se[i] = d->site_xpos[spair2[p_i]*3+i];  // position of the site in the environment
          n_cs[i] = csn[p_i*3+i];                   // contact surface normal
        }
        v_re  = p_se - p_sr;                        // vector from the end effector to the environment
        // distance
        phi_e = v_re.norm();                        // Euclidean distance between the end effector and the environment
        phi_n = v_re.dot(n_cs);                     // normal distance between the end effector and the environment
        zeta  = tanh(phi_r*phi_e);                    // semisphere based on the Euclidean distance
        phi_c = zeta*phi_e + (1-zeta)*phi_n;        // combined distance
        // normal force in the contact frame
        fn = kcon[p_i]*exp(-acon[p_i]*phi_c);
        // contact generalized in the world frame
        lambda = fn*n_cs;
        // loop for each free body
        for( int f_i=0; f_i<nfree; f_i++ )
        {
          for( int i=0; i<3; i++ )
          {
            p_fb[i] = d->qpos[jfree[f_i]+i];        // position of the center of mass of the free body
          }
          v_ef = p_cf - p_se;                       // vector from the end effector to the center of mass of free body
          // wrench on the free body due to the contact p_i: [lambda; cross(v_ef, lambda)]
          h[f_i*6+0] += lambda[0];
          h[f_i*6+1] += lambda[1];
          h[f_i*6+2] += lambda[2];
          h[f_i*6+3] += -v_ef[2]*lambda[1] + v_ef[1]*lambda[2];
          h[f_i*6+4] += v_ef[2]*lambda[0]  - v_ef[0]*lambda[2];
          h[f_i*6+5] += -v_ef[1]*lambda[0] + v_ef[0]*lambda[1];
        }
    }
    return h;
}
