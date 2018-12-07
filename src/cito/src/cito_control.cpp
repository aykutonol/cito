// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_CONTROl class consists of functions for calculating and setting control-
// and contact-related variables such as joint torques and external forces on
// the bodies.

// ***** CLASS TYPE ************************************************************
// Physics engine specific

#include "cito_control.h"

// ***** CONSTRUCTOR ***********************************************************


// ***** FUNCTIONS *************************************************************
// void CitoControl::takeFullStep(const mjModel* m, mjData* dmain, const ctrlVec_t ui, mjtNum* dXd)
// {
//     mjData* d;
//     d = mj_makeData(m);
//     // copy state and control from dmain to thread-specific d
//     d->time = dmain->time;
//     mju_copy(d->qpos, dmain->qpos, m->nq);
//     mju_copy(d->qvel, dmain->qvel, m->nv);
//     mju_copy(d->qacc, dmain->qacc, m->nv);
//     mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
//     mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
//     mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
//     mju_copy(d->ctrl, dmain->ctrl, m->nu);
//     // mju_copy(d->userdata, dmain->userdata, m->nuserdata);
//
//     // run full computation at center point (usually faster than copying dmain)
//     mj_forward(m, d);
//     sc.setControl(m, d, ui);
//
//     // take tc/dt steps
//     for( int j=0; j<ndpc; j++ )
//     {
//       // initialize the step
//       mj_step1(m, d);
//       // set ctrl and xfrc
//       sc.setControl(m, d, ui);
//       // complete the step
//       mj_step2(m, d);
//     }
//
//     // get new state
//     stateVec_t dXtemp; dXtemp.setZero();
//     dXtemp = this->getState(m, d);
//     mju_copy(dXd, dXtemp.data(), 2*m->nv);
//
//     // delete data
//     mj_deleteData(d);
// }

//setControl: sets generalized forces on joints and free bodies
void CitoControl::setControl(mjData* d, const ctrlVec_t u)
{
    // set control given the control input
    for( int i=0; i<params::nact; i++ )
    {
      d->ctrl[i] = u[i] + d->qfrc_bias[params::jact[i]];
    }
    // contact model
    hcon = this->contactModel(d, u);
    // set external forces on the free bodies
    for( int i=0; i<params::nfree; i++ )
    {
      for( int j=0; j<6; j++ )
      {
        d->xfrc_applied[params::bfree[i]*6+j] = hcon[i*6+j];
      }
    }
}

// contactModel: returns contact wrench given kcon
Eigen::Matrix<double, 6*params::nfree, 1> CitoControl::contactModel(const mjData* d, const ctrlVec_t u)
{
    h.setZero();
    // loop for each contact pair
    for( int p_i=0; p_i<params::npair; p_i++ )
    {
        // vectors in the world frame
        for( int i=0; i<3; i++ )
        {
          p_sr[i] = d->site_xpos[params::spair1[p_i]*3+i];  // position of the site on the robot
          p_se[i] = d->site_xpos[params::spair2[p_i]*3+i];  // position of the site in the environment
          n_cs[i] = params::csn[p_i*3+i];                   // contact surface normal
        }
        v_re  = p_se - p_sr;                                // vector from the end effector to the environment
        // distance
        phi_e = v_re.norm();                                // Euclidean distance between the end effector and the environment
        phi_n = v_re.dot(n_cs);                             // normal distance between the end effector and the environment
        zeta  = tanh(params::phi_r*phi_e);                  // semisphere based on the Euclidean distance
        phi_c = zeta*phi_e + (1-zeta)*phi_n;                // combined distance
        // normal force in the contact frame
        fn = u[params::nact+p_i]*exp(-params::acon[p_i]*phi_c);
        // contact generalized in the world frame
        lambda = fn*n_cs;
        // loop for each free body
        for( int f_i=0; f_i<params::nfree; f_i++ )
        {
          for( int i=0; i<3; i++ )
          {
            p_bf[i] = d->qpos[params::jfree[f_i]+i];        // position of the center of mass of the free body
          }
          v_ef = p_bf - p_se;                               // vector from the end effector to the free body
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
