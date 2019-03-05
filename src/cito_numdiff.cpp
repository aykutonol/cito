#include "cito_numdiff.h"

// ***** DESCRIPTION ***********************************************************
// CitoNumDiff class defines functions for numerical differentiation of MuJoCo
// dynamics including the forces imposed by the contact model.

// ***** CONSTRUCTOR ***********************************************************
CitoNumDiff::CitoNumDiff(const mjModel* model) : m(model), cc(model) {}
// ***** FUNCTIONS *************************************************************
// copyTakeStep: sets xNew to the integration of data given a control input
void CitoNumDiff::copyTakeStep(const mjData* dMain, const ctrlVec u, mjtNum* xNew)
{
    // create new data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, NV);
    mju_copy(d->qacc, dMain->qacc, NV);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, NV);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, NV);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dMain)
    mj_forward(m, d);
    cc.setControl(d, u);
    // take a full control step (i.e., tc/dt steps)
    cc.takeStep(d, u, false);
    // get new state
    xNewTemp.setZero();
    xNewTemp = cc.getState(d);
    mju_copy(xNew, xNewTemp.data(), N);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite-difference computation
void CitoNumDiff::hardWorker(const mjData* dMain, const ctrlVec uMain, mjtNum* deriv)
{
    // create data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, NV);
    mju_copy(d->qacc, dMain->qacc, NV);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, NV);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, NV);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
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
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data());
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
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
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data());
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
        // compute column i of dx/dqpos
        for( int j=0; j<N; j++ )
        {
            deriv[i*N + j] = (xNewP[j] - xNewN[j])/(2*eps);
        }
    }
    // finite-difference over velocities
    for( int i=0; i<NV; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data());
        // perturb velocity
        d->qvel[i] = dMain->qvel[i]-eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data());
        // undo perturbation
        d->qvel[i] = dMain->qvel[i];
        // compute column i of dx/dqvel
        for( int j=0; j<N; j++ )
        {
            deriv[N*NV + i*N + j] = (xNewP[j] - xNewN[j])/(2*eps);
        }
    }
    // finite-difference over control variables
    // copy uMain to utemp for perturbations
    utemp = uMain;
    for( int i=0; i<M; i++ )
    {
        // perturbation in the positive direction
        utemp[i] += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, utemp, xNewP.data());
        // perturbation in the negative direction
        utemp[i] -= 2*eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, utemp, xNewN.data());
        // compute column i of dx/du
        for( int j=0; j<N; j++ )
        {
            deriv[N*N + i*N + j] = (xNewP[j] - xNewN[j])/(2*eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

// linDyn: calculates derivatives of the state and control trajectories
void CitoNumDiff::linDyn(const mjData* dMain, const ctrlVec uMain, mjtNum* Fxd, mjtNum* Fud)
{
    mjtNum* deriv = (mjtNum*) mju_malloc(sizeof(mjtNum)*N*(N+M));
    this->hardWorker( dMain, uMain, deriv);
    mju_copy(Fxd, deriv, N*N);
    mju_copy(Fud, deriv+N*N, N*M);
    mju_free(deriv);
}
