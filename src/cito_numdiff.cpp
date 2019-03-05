#include "cito_numdiff.h"

// ***** DESCRIPTION ***********************************************************
// CitoNumDiff class defines functions for numerical differentiation of MuJoCo
// dynamics including the forces imposed by the contact model.

// ***** CONSTRUCTOR ***********************************************************
CitoNumDiff::CitoNumDiff(const mjModel* model) : m(model), cc(model)
{

}

// ***** FUNCTIONS *************************************************************
// copyTakeStep: sets newX to the integration of data given a control input
void CitoNumDiff::copyTakeStep(const mjData* dmain, const ctrlVec_t u, mjtNum* newX)
{
    // create new data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dmain to d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, NV);
    mju_copy(d->qacc, dmain->qacc, NV);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, NV);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, NV);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dmain)
    mj_forward(m, d);
    cc.setControl(d, u);
    // take a full control step (i.e., tc/dt steps)
    cc.takeStep(d, u, false);
    // get new state
    newXtemp.setZero();
    newXtemp = cc.getState(d);
    mju_copy(newX, newXtemp.data(), N);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite-difference computation
void CitoNumDiff::hardWorker(const mjData* dmain, const ctrlVec_t umain, mjtNum* deriv)
{
    // create data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dmain to d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, NV);
    mju_copy(d->qacc, dmain->qacc, NV);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, NV);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, NV);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);
    // finite-difference over positions
    for( int i=0; i<NV; i++ )
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
        this->copyTakeStep(d, umain, newXp.data());
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
        this->copyTakeStep(d, umain, newXn.data());
        // undo perturbation
        mju_copy(d->qpos, dmain->qpos, m->nq);
        // compute column i of dx/dqpos
        for( int j=0; j<N; j++ )
        {
            deriv[i*N + j] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // finite-difference over velocities
    for( int i=0; i<NV; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        newXp.setZero();
        this->copyTakeStep(d, umain, newXp.data());
        // perturb velocity
        d->qvel[i] = dmain->qvel[i]-eps;
        // get the negative perturbed state
        newXn.setZero();
        this->copyTakeStep(d, umain, newXn.data());
        // undo perturbation
        d->qvel[i] = dmain->qvel[i];
        // compute column i of dx/dqvel
        for( int j=0; j<N; j++ )
        {
            deriv[N*NV + i*N + j] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // finite-difference over control variables
    // copy umain to utemp for perturbations
    utemp = umain;
    for( int i=0; i<M; i++ )
    {
        // perturbation in the positive direction
        utemp[i] += eps;
        // get the positive perturbed state
        newXp.setZero();
        this->copyTakeStep(d, utemp, newXp.data());
        // perturbation in the negative direction
        utemp[i] -= 2*eps;
        // get the negative perturbed state
        newXn.setZero();
        this->copyTakeStep(d, utemp, newXn.data());
        // compute column i of dx/du
        for( int j=0; j<N; j++ )
        {
            deriv[N*N + i*N + j] = (newXp[j] - newXn[j])/(2*eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

// linDyn: calculates sensitivies of the states to changes in states and controls
void CitoNumDiff::linDyn(const mjData* dmain, const ctrlVec_t umain, mjtNum* Fxd, mjtNum* Fud)
{
    mjtNum* deriv = (mjtNum*) mju_malloc(sizeof(mjtNum)*N*(N+M));
    this->hardWorker(dmain, umain, deriv);
    mju_copy(Fxd, deriv, N*N);
    mju_copy(Fud, deriv+N*N, N*M);
    mju_free(deriv);
}
