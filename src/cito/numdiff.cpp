// ***** DESCRIPTION ***********************************************************
// NumDiff defines methods for numerical differentiation of MuJoCo
// dynamics including the forces imposed by the contact model.

#include "cito/numdiff.h"

// ***** CONSTRUCTOR ***********************************************************
NumDiff::NumDiff(const mjModel *m_, Params *cp_, Control *cc_) : m(m_), cp(cp_), cc(cc_)
{
    // initialize Eigen variables
    xNewTemp.resize(cp->n);
    xNewP.resize(cp->n);
    xNewN.resize(cp->n);
    uTemp.resize(cp->m);
}
// ***** FUNCTIONS *************************************************************
// copyTakeStep: sets xNew to the integration of data given a control input
void NumDiff::copyTakeStep(const mjData *dMain, const eigVd &u, double *xNew, double compensateBias)
{
    // create new data
    mjData *d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6 * m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // run full computation at center point (usually faster than copying dMain)
    mj_forward(m, d);
    cc->setControl(d, u, compensateBias);
    // take a full control step (i.e., tc/dt steps)
    cc->takeStep(d, u, 0, compensateBias);
    // get new state
    xNewTemp.setZero();
    xNewTemp = cp->getState(d);
    mju_copy(xNew, xNewTemp.data(), cp->n);
    // delete data
    mj_deleteData(d);
}

// hardWorker: for full, slow finite-difference computation
void NumDiff::hardWorker(const mjData *dMain, const eigVd &uMain, double *deriv, double compensateBias)
{
    // create data
    mjData *d;
    d = mj_makeData(m);
    // copy state and control from dMain to d
    d->time = dMain->time;
    mju_copy(d->qpos, dMain->qpos, m->nq);
    mju_copy(d->qvel, dMain->qvel, m->nv);
    mju_copy(d->qacc, dMain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dMain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dMain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dMain->xfrc_applied, 6 * m->nbody);
    mju_copy(d->ctrl, dMain->ctrl, m->nu);
    // finite-difference over positions
    for (int i = 0; i < m->nv; i++)
    {
        // get joint id for this dof
        int jID = m->dof_jntid[i];
        // apply quaternion or simple perturbation
        if (cp->quatAdr[i] >= 0)
        {
            mjtNum angvel[3] = {0, 0, 0};
            angvel[cp->dofAdr[i]] = eps;
            mju_quatIntegrate(d->qpos + cp->quatAdr[i], angvel, 1);
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
        if (cp->quatAdr[i] >= 0)
        {
            mjtNum angvel[3] = {0, 0, 0};
            angvel[cp->dofAdr[i]] = -eps;
            mju_quatIntegrate(d->qpos + cp->quatAdr[i], angvel, 1);
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
        for (int j = 0; j < cp->n; j++)
        {
            deriv[i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // finite-difference over velocities
    for (int i = 0; i < m->nv; i++)
    {
        // perturb velocity
        d->qvel[i] += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uMain, xNewP.data(), compensateBias);
        // perturb velocity
        d->qvel[i] = dMain->qvel[i] - eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uMain, xNewN.data(), compensateBias);
        // undo perturbation
        d->qvel[i] = dMain->qvel[i];
        // compute column i of dx/dqvel
        for (int j = 0; j < cp->n; j++)
        {
            deriv[cp->n * m->nv + i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // finite-difference over control variables
    // copy uMain to uTemp for perturbations
    uTemp = uMain;
    for (int i = 0; i < cp->m; i++)
    {
        // perturbation in the positive direction
        uTemp(i) += eps;
        // get the positive perturbed state
        xNewP.setZero();
        this->copyTakeStep(d, uTemp, xNewP.data(), compensateBias);
        // perturbation in the negative direction
        uTemp(i) -= 2 * eps;
        // get the negative perturbed state
        xNewN.setZero();
        this->copyTakeStep(d, uTemp, xNewN.data(), compensateBias);
        // compute column i of dx/du
        for (int j = 0; j < cp->n; j++)
        {
            deriv[cp->n * cp->n + i * cp->n + j] = (xNewP(j) - xNewN(j)) / (2 * eps);
        }
    }
    // delete data
    mj_deleteData(d);
}

// linDyn: calculates derivatives of the state and control trajectories
void NumDiff::linDyn(const mjData *dMain, const eigVd &uMain, double *Fxd, double *Fud, double compensateBias)
{
    // TODO: consider doing the memory allocation/freeing in the constructor/destructor
    double *deriv = (double *)mju_malloc(sizeof(double) * cp->n * (cp->n + cp->m));
    this->hardWorker(dMain, uMain, deriv, compensateBias);
    mju_copy(Fxd, deriv, cp->n * cp->n);
    mju_copy(Fud, deriv + cp->n * cp->n, cp->n * cp->m);
    mju_free(deriv);
}
