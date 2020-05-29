#include "cito_numdiff.h"

// ***** DESCRIPTION ***********************************************************
// CitoNumDiff class defines functions for numerical differentiation of MuJoCo
// dynamics including the forces imposed by the contact model.

// ***** CONSTRUCTOR ***********************************************************
CitoNumDiff::CitoNumDiff(const mjModel* model) : m(model), cp(model), cc(model)
{
    // initialize Eigen variables
    xNewTemp.resize(cp.n); xNewP.resize(cp.n); xNewN.resize(cp.n);
    uTemp.resize(cp.m);
}
// ***** FUNCTIONS *************************************************************
// copyTakeStep: sets xNew to the integration of data given a control input
void CitoNumDiff::copyTakeStep(const mjData* dMain, const eigVd u, mjtNum* xNew, double compensateBias)
{
    // create new data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dMain)
    mj_forward(m, d);
    cc.setControl(d, u, compensateBias);
    // take a full control step (i.e., tc/dt steps)
    cc.takeStep(d, u, false, compensateBias);
    // get new state
    xNewTemp.setZero();
    xNewTemp = cc.getState(d);
    mju_copy(xNew, xNewTemp.data(), cp.n);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite-difference computation
void CitoNumDiff::hardWorker(const mjData* dMain, const eigVd uMain, mjtNum* deriv, double compensateBias)
{
    // create data
    mjData* d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // finite-difference over positions
    for( int i=0; i<m->nv; i++ )
    {
        // get joint id for this dof
        int jID = m->dof_jntid[i];
        // apply quaternion or simple perturbation
        if( cp.quatAdr[i]>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[cp.dofAdr[i]] = eps;
            mju_quatIntegrate(d->qpos+cp.quatAdr[i], angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jID] + i - m->jnt_dofadr[jID]] += eps;
        }
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data(), compensateBias);
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
        // apply quaternion or simple perturbation
        if( cp.quatAdr[i]>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[cp.dofAdr[i]] = -eps;
            mju_quatIntegrate(d->qpos+cp.quatAdr[i], angvel, 1);
        }
        else
        {
            d->qpos[m->jnt_qposadr[jID] + i - m->jnt_dofadr[jID]] -= eps;
        }
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data(), compensateBias);
        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);
        // compute column i of dx/dqpos
        for( int j=0; j<cp.n; j++ )
        {
            deriv[i*cp.n + j] = (xNewP(j) - xNewN(j))/(2*eps);
        }
    }
    // finite-difference over velocities
    for( int i=0; i<m->nv; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data(), compensateBias);
        // perturb velocity
        d->qvel[i] = dMain->qvel[i]-eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data(), compensateBias);
        // undo perturbation
        d->qvel[i] = dMain->qvel[i];
        // compute column i of dx/dqvel
        for( int j=0; j<cp.n; j++ )
        {
            deriv[cp.n*m->nv + i*cp.n + j] = (xNewP(j) - xNewN(j))/(2*eps);
        }
    }
    // finite-difference over control variables
    // copy uMain to uTemp for perturbations
    uTemp = uMain;
    for( int i=0; i<cp.m; i++ )
    {
        // perturbation in the positive direction
        uTemp(i) += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uTemp, xNewP.data(), compensateBias);
        // perturbation in the negative direction
        uTemp(i) -= 2*eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uTemp, xNewN.data(), compensateBias);
        // compute column i of dx/du
        for( int j=0; j<cp.n; j++ )
        {
            deriv[cp.n*cp.n + i*cp.n + j] = (xNewP(j) - xNewN(j))/(2*eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

// linDyn: calculates derivatives of the state and control trajectories
void CitoNumDiff::linDyn(const mjData* dMain, const eigVd uMain, mjtNum* Fxd, mjtNum* Fud, double compensateBias)
{
    mjtNum* deriv = (mjtNum*) mju_malloc(sizeof(mjtNum)*cp.n*(cp.n+cp.m));
    this->hardWorker( dMain, uMain, deriv, compensateBias);
    mju_copy(Fxd, deriv, cp.n*cp.n);
    mju_copy(Fud, deriv+cp.n*cp.n, cp.n*cp.m);
    mju_free(deriv);
}

// worker: performs fast finite-difference computation
void CitoNumDiff::worker(const mjData* dMain, mjtNum* deriv)
{
    // create data
    mjData* d;
    d = mj_makeData(m);
    // allocate stack space for result at center
    mjMARKSTACK
    mjtNum* center = mj_stackAlloc(d, m->nv);
    mjtNum* warmstart = mj_stackAlloc(d, m->nv);

    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);

    // run full computation at center point (usually faster than copying dMain)
    mj_forward(m, d);
    // extra solver iterations to improve warmstart (qacc) at center point
    for( int rep=1; rep<this->nwarmup; rep++ )
        mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

    // select output from forward or inverse dynamics
    mjtNum* output = d->qacc;

    // save output for center point and warmstart (needed in forward only)
    mju_copy(center, output, m->nv);
    mju_copy(warmstart, d->qacc_warmstart, m->nv);

    // select target vector and original vector for force or acceleration derivative
    mjtNum* target = d->qfrc_applied;
    const mjtNum* original = dMain->qfrc_applied;

    // finite-difference over force or acceleration: skip = mjSTAGE_VEL
    for( int i=0; i<m->nv; i++ )
    {
        // perturb selected target
        target[i] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

        // undo perturbation
        target[i] = original[i];

        // compute column i of derivative 2
        for( int j=0; j<m->nv; j++ )
            deriv[2*m->nv*m->nv + i + j*m->nv] = (output[j] - center[j])/eps;
    }

    // finite-difference over velocity: skip = mjSTAGE_POS
    for( int i=0; i<m->nv; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_POS, 1);

        // undo perturbation
        d->qvel[i] = dMain->qvel[i];

        // compute column i of derivative 1
        for( int j=0; j<m->nv; j++ )
            deriv[m->nv*m->nv + i + j*m->nv] = (output[j] - center[j])/eps;
    }

    // finite-difference over position: skip = mjSTAGE_NONE
    for( int i=0; i<m->nv; i++ )
    {
        // get joint id for this dof
        int jid = m->dof_jntid[i];

        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        int quatadr = -1, dofpos = 0;
        if( m->jnt_type[jid]==mjJNT_BALL )
        {
            quatadr = m->jnt_qposadr[jid];
            dofpos = i - m->jnt_dofadr[jid];
        }
        else if( m->jnt_type[jid]==mjJNT_FREE && i>=m->jnt_dofadr[jid]+3 )
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
            d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_NONE, 1);

        // undo perturbation
        mju_copy(d->qpos, dMain->qpos, m->nq);

        // compute column i of derivative 0
        for( int j=0; j<m->nv; j++ )
            deriv[i + j*m->nv] = (output[j] - center[j])/eps;
    }
    // delete data
    mj_deleteData(d);
}