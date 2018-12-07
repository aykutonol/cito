// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

// ***** DESCRIPTION ***********************************************************
// CITO_NUMDIFF class consists of functions for numerical differentiation of the
// MuJoCo dynamics including the contact forces imposed by the contact model.

// ***** CLASS TYPE ************************************************************
// Physics engine specific

#include "cito_numdiff.h"

// ***** CONSTRUCTOR ***********************************************************
CitoNumDiff::CitoNumDiff(const mjModel* model, const mjData* data, const ctrlVec_t control)
: m(model), dmain(data), umain(control)
{
}

// ***** FUNCTIONS *************************************************************
void CitoNumDiff::copyTakeStep(mjtNum* newX, ctrlVec_t u)
{
    // create new data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dmain to d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, m->nv);
    mju_copy(d->qacc, dmain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dmain)
    mj_forward(m, d);
    cc.setControl(m, d, umain);
    // take a full control step (i.e., tc/dt steps)
    for( int j=0; j<ndpc; j++ )
    {
      mj_step1(m, d);
      cc.setControl(d, u);
      mj_step2(m, d);
    }
    // get new state
    newXtemp.setZero();
    newXtemp = this->getState(d);
    mju_copy(newX, newXtemp.data(), 2*m->nv);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite difference computation
void CitoNumDiff::hardWorker(mjtNum* deriv)
{
    // make data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dmain to d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, m->nv);
    mju_copy(d->qacc, dmain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);
    // copy umain to up for perturbations
    up = umain;
    // finite-difference over control variables
    for( int i=0; i<M; i++ )
    {
        // perturbation in the positive direction
        up[i] += eps;
        // get the positive perturbed state
        newXp.setZero();
        this->copyTakeStep(newXp.data(), up);
        // undo perturbation
        up[i] = umain[i];
        // perturbation in the negative direction
        up[i] -= eps;
        // get the negative perturbed state
        newXn.setZero();
        this->copyTakeStep(newXn.data(), up);
        // undo perturbation
        up[i] = umain[i];
        // compute column i of dx/du
        for( int j=0; j<2*m->nv; j++ )
        {
          deriv[4*m->nv*m->nv + i + j*M] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // finite-difference over velocities
    for( int i=0; i<m->nv; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        newXp.setZero();
        this->copyTakeStep(newXp.data(), umain);
        // undo perturbation
        d->qvel[i] = dmain->qvel[i];
        // perturb velocity
        d->qvel[i] -= eps;
        // get the negative perturbed state
        newXn.setZero();
        this->copyTakeStep(newXn.data(), umain);
        // undo perturbation
        d->qvel[i] = dmain->qvel[i];
        // compute column i of dx/dqvel
        for( int j=0; j<2*m->nv; j++ )
        {
          deriv[2*m->nv*m->nv + i + j*m->nv] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // finite-difference over positions
    for( int i=0; i<m->nv; i++ )
    {
        // get joint id for this dof
        int jid = m->dof_jntid[i];
        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        int quatadr = -1, dofpos = 0;
        if( m->jnt_type[jid]==mjJNT_FREE && i>=m->jnt_dofadr[jid]+3 )
        {
            quatadr = m->jnt_qposadr[jid] + 3;
            dofpos = i - m->jnt_dofadr[jid] - 3;
        }
        // apply quaternion or simple perturbation
        if( quatadr>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[dofpos] = eps;
            mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] += eps;
        }
        // get the positive perturbed state
        newXp.setZero();
        this->copyTakeStep(newXp.data(), umain);
        // undo perturbation
        mju_copy(d->qpos, dmain->qpos, m->nq);
        // apply quaternion or simple perturbation
        if( quatadr>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[dofpos] = -eps;
            mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] -= eps;
        }
        // get the negative perturbed state
        newXn.setZero();
        this->copyTakeStep(newXn.data(), umain);
        // undo perturbation
        mju_copy(d->qpos, dmain->qpos, m->nq);
        // compute column i of dx/dqpos
        for( int j=0; j<2*m->nv; j++ )
        {
            deriv[i + j*m->nv] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

// getState function converts free joints' quaternions to Euler angles so that
// the dimensionality of the state vector is 2*nv instead of nq+nv
stateVec_t CitoNumDiff::getState(const mjData* d)
{
    stateVec_t x; x.setZero();
    int free_count = 0;
    if ( m->nq != m->nv )
    {
        for ( int i=0; i<m->nq; i++ )
        {
            int jid = m->dof_jntid[i];
            if( m->jnt_type[jid]==mjJNT_FREE )
            {
                mju_copy(x.block<3,1>(i,0).data(), d->qpos + i, 3);
                Eigen::Matrix<mjtNum, 4, 1> obj_q;
                mju_copy(obj_q.data(), d->qpos + i + 3, 4);
                // calculate euler angles from the quaternion
                x(i+3) = atan2(2*(obj_q[0]*obj_q[1]+obj_q[2]*obj_q[3]), 1-2*(pow(obj_q[1],2)+pow(obj_q[2],2)));
                x(i+4) = asin(2*(obj_q[0]*obj_q[2]-obj_q[3]*obj_q[1]));
                x(i+5) = atan2(2*(obj_q[0]*obj_q[3]+obj_q[1]*obj_q[2]), 1-2*(pow(obj_q[2],2)+pow(obj_q[3],2)));
                i += 6;             // skip next 6 position indices
                free_count += 1;    // free joint counter
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
    mju_copy(x.data() + m->nv, d->qvel, m->nv);

    return x;
}


// linDyn function calculates sensitivies of the states to changes in states and controls
void CitoNumDiff::linDyn(mjtNum* Fxd, mjtNum* Fud)
{
    mjtNum* deriv = (mjtNum*) mju_malloc(6*sizeof(mjtNum)*m->nv*m->nv);
    this->hardWorker(deriv);
    mju_copy(Fxd, deriv, 4*m->nv*m->nv);
    mju_copy(Fud, deriv+4*m->nv*m->nv, 2*m->nv*M);
    mju_free(deriv);
}
