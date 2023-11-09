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
//    std::cout << "model nv: " << m->nv << std::endl;
//    std::cout << "model nq: " << m->nq << std::endl;
//    std::cout << "model nu: " << m->nu << std::endl;
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

std::vector<int> NumDiff::generateKeypoints(derivative_interpolator di, const eigMd X, int horizon){
    std::vector<int> keypoints;
    keypoints.push_back(0);

    if(di.keyPoint_method == "set_interval"){
        for(int i = 1; i < horizon; i++){
            if(i % di.min_n == 0){
                keypoints.push_back(i);
            }
        }

    }
    else if(di.keyPoint_method == "adaptive_jerk"){

    }
    else if(di.keyPoint_method == "magvel_change"){

    }
    else if(di.keyPoint_method == "iterative_error"){

    }
    else{
        std::cout << "keyPoint_method not recognized" << std::endl;
    }

    // Check if last keypoint is the horizon
    if(keypoints.back() != horizon){
        keypoints.push_back(horizon);
    }

    return keypoints;
}

void NumDiff::interpolateDerivs(std::vector<int> keypoints, eigTd &Fxd, eigTd &Fud, int horizon){
    int num_keypoints = keypoints.size();
    int dof = Fxd[0].rows() / 2;
    int num_ctrl = Fud[0].cols();

    // Loop over the keypoints
    for(int i = 0; i < num_keypoints - 1; i++){
        int start_index = keypoints[i];
        int end_index = keypoints[i+1];
        int interval = end_index - start_index;

        eigMd Fxd_start = Fxd[start_index];
        eigMd Fxd_end = Fxd[end_index];
        eigMd Fxd_add = (Fxd_end - Fxd_start) / interval;

        eigMd Fud_start = Fud[start_index];
        eigMd Fud_end = Fud[end_index];
        eigMd Fud_add = (Fud_end - Fud_start) / interval;

        for(int j = 1; j < interval; j++){
            Fxd[start_index + j] = Fxd_start + Fxd_add * j;
            Fud[start_index + j] = Fud_start + Fud_add * j;
        }
    }

//    // Interpolate Fud
//    for(int i = 0; i < num_keypoints - 1; i++){
//        int start = keypoints[i];
//        int end = keypoints[i+1];
//        int num_intervals = end - start;
//        for(int j = 0; j < num_intervals; j++){
//            double alpha = (double)j / (double)num_intervals;
//            Fud[start + j] = (1 - alpha) * Fud[start] + alpha * Fud[end];
//        }
//    }
//
//    // Fill in the rest of the trajectory
//    for(int i = 0; i < horizon; i++){
//        if(Fxd[i].isZero()){
//            Fxd[i] = Fxd[i-1];
//        }
//        if(Fud[i].isZero()){
//            Fud[i] = Fud[i-1];
//        }
//    }
}

void NumDiff::saveLinearisation(const std::string file_prefix, eigTd Fxd, eigTd Fud, int horizon){
    std::string projectParentPath = __FILE__;
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    projectParentPath = projectParentPath.substr(0, projectParentPath.find_last_of("/\\"));
    std::string rootPath = projectParentPath + "/savedTrajecInfo/" + file_prefix;

    std::string filename = rootPath + "_A_matrices.csv";
    std::cout << "filename: " << filename << std::endl;
    std::ofstream fileOutput;
    fileOutput.open(filename);

    int dof = Fxd[0].rows() / 2;
    int num_ctrl = Fud[0].cols();

    // trajectory length
    for(int i = 0; i < horizon - 1; i++){
        // Row
        for(int j = 0; j < (dof); j++){
            // Column
            for(int k = 0; k < (2 * dof); k++){
//                std::cout << Fxd[i](j, k) << ",";
                fileOutput << Fxd[i](j + dof, k) << ",";
            }

        }
        fileOutput << std::endl;
    }

    fileOutput.close();
}
