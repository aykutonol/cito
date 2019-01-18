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
CitoControl::CitoControl(const mjModel* model) : m(model) {}

CitoControl::~CitoControl()
{
    delete []qpos_lb;       delete []qpos_ub;
    delete []tau_lb;        delete []tau_ub;
    delete []isJFree;       delete []isAFree;
}

// ***** FUNCTIONS *************************************************************
// takeStep: takes a full control step given a control input
void CitoControl::takeStep(mjData*d, const ctrlVec_t u)
{
    for( int i=0; i<params::ndpc; i++ )
    {
        mj_step1(m, d);
        this->setControl(d, u);
        mj_step2(m, d);
    }
}

// setControl: sets generalized forces on joints and free bodies
void CitoControl::setControl(mjData* d, const ctrlVec_t u)
{
    // set control given the control input
    for( int i=0; i<NU; i++ )
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

// contactModel: returns contact wrench given finalPose
Eigen::Matrix<double, 6*params::nfree, 1> CitoControl::contactModel(const mjData* d, const ctrlVec_t u)
{
    h.setZero();
    // loop for each contact pair
    for( int p_i=0; p_i<NPAIR; p_i++ )
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
        fn = u[NU+p_i]*exp(-params::acon[p_i]*phi_c);
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

// getState function converts free joints' quaternions to Euler angles so that
// the dimensionality of the state vector is 2*nv instead of nq+nv
stateVec_t CitoControl::getState(const mjData* d)
{
    x.setZero();
    int free_count = 0;
    if ( m->nq != NV )
    {
        for ( int i=0; i<m->nq; i++ )
        {
            int jid = m->dof_jntid[i];
            if( m->jnt_type[jid]==mjJNT_FREE )
            {
                mju_copy(x.block<3,1>(i,0).data(), d->qpos+i, 3);
                mju_copy(jfree_quat.data(), d->qpos+i+3, 4);
                // calculate the Euler angles from the quaternion
                x(i+3) = atan2(2*(jfree_quat[0]*jfree_quat[1]+jfree_quat[2]*jfree_quat[3]), 1-2*(pow(jfree_quat[1],2)+pow(jfree_quat[2],2)));
                x(i+4) =  asin(2*(jfree_quat[0]*jfree_quat[2]-jfree_quat[3]*jfree_quat[1]));
                x(i+5) = atan2(2*(jfree_quat[0]*jfree_quat[3]+jfree_quat[1]*jfree_quat[2]), 1-2*(pow(jfree_quat[2],2)+pow(jfree_quat[3],2)));
                i += 6;             // proceed to next joint
                free_count++;       // free joint counter
            }
            else
            {
                x[i-free_count] = d->qpos[i];
            }
        }
    }
    else
    {
        mju_copy(x.data(), d->qpos, m->nq);
    }
    // get the velocities
    mju_copy(x.data() + NV, d->qvel, NV);

    return x;
}

// getBounds: gets bounds on joint positions, actuator forces from the model
void CitoControl::getBounds()
{
    for( int i=0; i<NV; i++ )
    {
        int jid =  m->dof_jntid[i];
        if( m->jnt_limited[jid] )
        {
            isJFree[i] = 0;
            qpos_lb[i] = m->jnt_range[jid*2];
            qpos_ub[i] = m->jnt_range[jid*2+1];
        }
        else
        {
            isJFree[i] = 1;
            qpos_lb[i] = 0;  // to be replaced by -infBnd in the initial guess
            qpos_ub[i] = 0;  // to be replaced by +infBnd in the initial guess
        }
    }
    for( int i=0; i<NU; i++ )
    {
        if( m->actuator_ctrllimited[i] )
        {
            isAFree[i] = 0;
            tau_lb[i]  = m->actuator_ctrlrange[i*2];
            tau_ub[i]  = m->actuator_ctrlrange[i*2+1];
        }
        else
        {
            isAFree[i] = 1;
            tau_lb[i]  = 0; // to be replaced by -infBnd in the initial guess
            tau_ub[i]  = 0; // to be replaced by +infBnd in the initial guess
        }
    }
}