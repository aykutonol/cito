// ***** DESCRIPTION ***********************************************************
// Control defines functions for calculating and setting control- and contact-
// related variables, e.g., joint torques and external forces on the free
// bodies.ss
#include "cito/control.h"

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
Control::Control(const mjModel *m_, Params *cp_) : m(m_), cp(cp_), sl(m_, cp_)
{
    // read contact model parameters
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/vscm.yaml");
    alpha = vscm["alpha"].as<double>();
    // bound variables
    qposLB = new double[m->nv];
    qposUB = new double[m->nv];
    tauLB = new double[m->nu];
    tauUB = new double[m->nu];
    isJFree = new int[m->nv];
    isAFree = new int[m->nu];
    // initialize Eigen variables
    h.resize(6 * cp->nFree);
    hCon.resize(6 * cp->nFree);
    // set FCL distance request options
    distReq.enable_nearest_points = true;
    distReq.enable_signed_distance = true;
    // create collision geometries/objects
    fcl::Transform3d tf0;
    tf0.translation().setZero();
    tf0.linear().setIdentity();
    for (int i = 0; i < cp->nPair; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if (collObjs[cp->sites[i][j]] == NULL)
                collObjs[cp->sites[i][j]] = new fcl::CollisionObjectd(createCollGeom(m, cp->sites[i][j]), tf0);
        }
    }
}
Control::~Control()
{
    // delete bound variables
    delete[] qposLB;
    delete[] qposUB;
    delete[] tauLB;
    delete[] tauUB;
    delete[] isJFree;
    delete[] isAFree;
    // delete collision objects
    for (auto it : collObjs)
        delete it.second;
}

// ***** FUNCTIONS *************************************************************
// getSiteTransform: returns the site pose from MuJoCo data
fcl::Transform3d Control::getSiteTransform(const mjData *d, int site_id)
{
    fcl::Transform3d tf;
    // get the position
    mju_copy3(tf.translation().data(), d->site_xpos + 3 * site_id);
    // get the rotation
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R;
    mju_copy(R.data(), d->site_xmat + 9 * site_id, 9);
    tf.linear() = R;
    return tf;
}

// createCollGeom creates an FCL collision object given a MuJoCo model, site no, and transform
std::shared_ptr<fcl::CollisionGeometryd> Control::createCollGeom(const mjModel *m, int site_id)
{
    std::shared_ptr<fcl::CollisionGeometryd> geom = NULL;
    if (m->site_type[site_id] == mjGEOM_SPHERE)
    {
        printf("\tCreating a spherical collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Sphered>(m->site_size[site_id * 3]);
    }
    else if (m->site_type[site_id] == mjGEOM_CYLINDER)
    {
        printf("\tCreating a cylindirical collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Cylinderd>(m->site_size[site_id * 3 + 1] * 2,
                                                m->site_size[site_id * 3 + 0]);
    }
    else if (m->site_type[site_id] == mjGEOM_BOX)
    {
        printf("\tCreating a prismatic collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Boxd>(m->site_size[site_id * 3 + 0] * 2,
                                           m->site_size[site_id * 3 + 1] * 2,
                                           m->site_size[site_id * 3 + 2] * 2);
    }
    else
    {
        printf("\t\033[0;31mCannot add site %d b/c geometry type %d is not implemented.\033[0m\n",
               site_id, m->site_type[site_id]);
    }
    return geom;
}

// calcDistance: returns FCL distance calculation result for each pair
std::vector<fcl::DistanceResultd> Control::calcDistance(const mjData *d)
{
    std::vector<fcl::DistanceResultd> distResults;
    // update collision objects' poses
    for (auto collObj : collObjs)
        collObj.second->setTransform(getSiteTransform(d, collObj.first));
    // FCL distance calculation
    for (int pair = 0; pair < cp->nPair; pair++)
    {
        distRes.clear();
        fcl::distance(collObjs[cp->sites[pair][0]], collObjs[cp->sites[pair][1]], distReq, distRes);
        distResults.push_back(distRes);
    }
    return distResults;
}

// takeStep: takes a full control step given a control input
void Control::takeStep(mjData *d, const eigVd &u, int save, double compensateBias)
{
    if (save > 0 && d->time <= 1e-6)
    {
        sl.writeData(d, save);
    }
    for (int i = 0; i < cp->ndpc; i++)
    {
        mj_step1(m, d);
        this->setControl(d, u, compensateBias);
        mj_step2(m, d);
        if (save > 0)
        {
            sl.writeData(d, save);
        }
    }
}

// setControl: sets generalized forces on joints and free bodies
void Control::setControl(mjData *d, const eigVd &u, double compensateBias)
{
    // set control given the control input
    for (int i = 0; i < m->nu; i++)
    {
        d->ctrl[i] = u(i) + compensateBias * d->qfrc_bias[cp->dAct[i]] - 0. * d->qfrc_constraint[cp->dAct[i]];
    }
    // contact model
    hCon.setZero();
    hCon = this->contactModel(d, u);
    // set external forces on the free bodies
    for (int i = 0; i < cp->nFree; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            d->xfrc_applied[cp->bFree[i] * 6 + j] = hCon(i * 6 + j);
        }
    }
}

// contactModel: returns contact wrench given current state and control input
eigVd Control::contactModel(const mjData *d, const eigVd &u)
{
    h.setZero();
    // calculate min. distance and nearest points for all contact pairs
    std::vector<fcl::DistanceResultd> distPairs = calcDistance(d);
    // loop for each contact pair
    for (int pair = 0; pair < cp->nPair; pair++)
    {
        // update contact surface normal
        mju_rotVecMat(nCS.data(), cp->unit_x, d->site_xmat + 9 * cp->sites[pair][1]);
        // calculate the normal force in the contact frame
        gamma = u(m->nu + pair) * exp(-alpha * distPairs[pair].min_distance);
        // contact the linear force in the world frame
        lambda = gamma * nCS;
        // map the force to the CoM of each free body
        for (int free_body = 0; free_body < cp->nFree; free_body++)
        {
            // get the free body's CoM position, i.e., joint position for a free joint
            mju_copy3(pCoM.data(), d->qpos + cp->pFree[free_body]);
            // calculate the vector from the CoM to the contact point in the environment
            // r = distPairs[pair].nearest_points[1] - pCoM;
            mju_copy3(sEnvPos.data(), d->site_xpos + 3 * cp->sites[pair][1]);
            r = sEnvPos - pCoM;
            // calculate the wrench at the CoM: [lambda; cross(vEF, lambda)]
            h.segment(free_body * 6, 3) += lambda;
            h.segment(free_body * 6 + 3, 3) += cp->skewCross(r, lambda);
        }
    }
    return h;
}

// getBounds: gets bounds on joint positions, actuator forces from the model
void Control::getBounds()
{
    for (int i = 0; i < m->nv; i++)
    {
        int jID = m->dof_jntid[i];
        if (m->jnt_limited[jID])
        {
            isJFree[i] = 0;
            qposLB[i] = m->jnt_range[jID * 2];
            qposUB[i] = m->jnt_range[jID * 2 + 1];
        }
        else
        {
            isJFree[i] = 1;
            qposLB[i] = 0; // to be replaced by -infBnd in the initial guess
            qposUB[i] = 0; // to be replaced by +infBnd in the initial guess
        }
    }
    for (int i = 0; i < m->nu; i++)
    {
        if (m->actuator_ctrllimited[i])
        {
            isAFree[i] = 0;
            tauLB[i] = m->actuator_ctrlrange[i * 2];
            tauUB[i] = m->actuator_ctrlrange[i * 2 + 1];
        }
        else
        {
            isAFree[i] = 1;
            tauLB[i] = 0; // to be replaced by -infBnd in the initial guess
            tauUB[i] = 0; // to be replaced by +infBnd in the initial guess
        }
    }
}