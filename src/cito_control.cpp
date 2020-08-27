// ***** DESCRIPTION ***********************************************************
// CitoControl class defines functions for calculating and setting control- and
// contact-related variables, i.e., joint torques and external forces on the free
// bodies.

#include "cito_control.h"


// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoControl::CitoControl(const mjModel* m_, CitoParams* cp_) : m(m_), cp(cp_), sl(m_, cp_)
{
    // read contact model parameters
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    alpha = vscm["alpha"].as<double>();
    // bound variables
    qposLB  = new double[m->nv];   qposUB  = new double[m->nv];
    tauLB   = new double[m->nu];   tauUB   = new double[m->nu];
    isJFree = new int[m->nv];      isAFree = new int[m->nu];
    // initialize Eigen variables
    h.resize(6*cp->nFree); hCon.resize(6*cp->nFree);
    x.resize(cp->n);
    // set FCL distance request options
    distReq.enable_nearest_points = true;
    distReq.enable_signed_distance = true;
    // create collision geometries/objects
    fcl::Transform3d tf0;
    tf0.translation().setZero(); tf0.linear().setIdentity();
    for( int i=0; i<cp->nPair; i++ )
    {
        for( int j=0; j<2; j++) {
            if(collObjs[cp->sites[i][j]]==NULL)
                collObjs[cp->sites[i][j]] = new fcl::CollisionObjectd(createCollGeom(m, cp->sites[i][j]), tf0);
        }
    }
}
CitoControl::~CitoControl()
{
    // delete bound variables
    delete[] qposLB;        delete[] qposUB;
    delete[] tauLB;         delete[] tauUB;
    delete[] isJFree;       delete[] isAFree;
    // delete collision objects
    for(auto it : collObjs)
        delete it.second;
}

// ***** FUNCTIONS *************************************************************
// getSiteTransform: returns the site pose from MuJoCo data
fcl::Transform3d CitoControl::getSiteTransform(const mjData* d, int site_id) {
    fcl::Transform3d tf;
    // get the position
    mju_copy3(tf.translation().data(), d->site_xpos+3*site_id);
    // get the rotation
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R;
    mju_copy(R.data(), d->site_xmat+9*site_id, 9);
    tf.linear() = R;
    return tf;
}

// createCollGeom creates an FCL collision object given a MuJoCo model, site no, and transform
std::shared_ptr<fcl::CollisionGeometryd> CitoControl::createCollGeom(const mjModel* m, int site_id) {
    std::shared_ptr<fcl::CollisionGeometryd> geom = NULL;
    if( m->site_type[site_id]==mjGEOM_SPHERE )
    {
        printf("\tCreating a spherical collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Sphered>(m->site_size[site_id*3]);
    }
    else if( m->site_type[site_id]==mjGEOM_CYLINDER )
    {
        printf("\tCreating a cylindirical collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Cylinderd>(m->site_size[site_id*3+1]*2,
                                                m->site_size[site_id*3+0]);
    }
    else if( m->site_type[site_id]==mjGEOM_BOX )
    {
        printf("\tCreating a prismatic collision object for site %d.\n", site_id);
        geom = std::make_shared<fcl::Boxd>(m->site_size[site_id*3+0]*2,
                                           m->site_size[site_id*3+1]*2,
                                           m->site_size[site_id*3+2]*2);
    }
    else
    {
        printf("\t\033[0;31mCannot add site %d b/c geometry type %d is not implemented.\033[0m\n",
               site_id, m->site_type[site_id]);
    }
    return geom;
}

// takeStep: takes a full control step given a control input
void CitoControl::takeStep(mjData* d, const eigVd u, bool save, double compensateBias)
{
    if( save ) { sl.writeData(d); }
    for( int i=0; i<cp->ndpc; i++ )
    {
        mj_step1(m, d);
        this->setControl(d, u, compensateBias);
        mj_step2(m, d);
        if( save ) { sl.writeData(d); }
    }
}

// setControl: sets generalized forces on joints and free bodies
void CitoControl::setControl(mjData* d, const eigVd u, double compensateBias)
{
    // set control given the control input
    for( int i=0; i<m->nu; i++ )
    {
      d->ctrl[i] = u(i) + compensateBias*d->qfrc_bias[cp->dAct[i]];
    }
    // contact model
    hCon.setZero();
    hCon = this->contactModel(d, u);
    // set external forces on the free bodies
    for( int i=0; i<cp->nFree; i++ )
    {
      for( int j=0; j<6; j++ )
      {
        d->xfrc_applied[cp->bFree[i]*6+j] = hCon(i*6+j);
      }
    }
}

// contactModel: returns contact wrench given current state and control input
eigVd CitoControl::contactModel(const mjData* d, const eigVd u)
{
    h.setZero();
    // update collision objects' poses
    for(auto it : collObjs)
        it.second->setTransform(getSiteTransform(d, it.first));
    // loop for each contact pair
    for( int pair=0; pair<cp->nPair; pair++ )
    {
        // contact surface normal
        mju_rotVecMat(nCS.data(), unit_x, d->site_xmat+9*cp->sites[pair][1]);
        // FCL distance calculation
        distRes.clear();
        fcl::distance(collObjs[cp->sites[pair][0]], collObjs[cp->sites[pair][1]], distReq, distRes);
        // normal force in the contact frame
        gamma = u(m->nu+pair)*exp(-alpha*distRes.min_distance);
        // contact generalized in the world frame
        lambda = gamma*nCS;
        // map the force to the CoM of each free body
        for( int free_body=0; free_body<cp->nFree; free_body++ )
        {
            // get the free body's CoM position, i.e., joint position for a free joint
            mju_copy3(pCoM.data(), d->qpos+cp->pFree[free_body]);
            // calculate the vector from the CoM to the contact point in the environment
            r = distRes.nearest_points[1] - pCoM;
            // calculate the wrench at the CoM: [lambda; cross(vEF, lambda)]
            h.segment(free_body*6, 3) += lambda;
            h.segment(free_body*6+3, 3) += cp->skewCross(r, lambda);
        }
    }
    return h;
}

// getState: converts free joints' quaternions to Euler angles so that
// the dimensionality of the state vector is 2*nv instead of nq+nv
eigVd CitoControl::getState(const mjData* d)
{
    x.setZero();
    int freeNo = 0;
    // get the positions
    for ( int i=0; i<m->nq; i++ )
    {
        int jID = m->dof_jntid[i-freeNo];
        if( m->jnt_type[jID] == mjJNT_FREE )
        {
            jFreeQuat.setZero();
            mju_copy3(x.segment(i, 3).data(), d->qpos+i);
            mju_copy4(jFreeQuat.data(), d->qpos+i+3);
            // convert quaternion into Euler angles
            x.segment(i+3, 3) = cp->quat2Euler(jFreeQuat);
            // count the free joint and proceed to next joint
            freeNo++;
            i += 6;
        }
        else
        {
            x(i-freeNo) = d->qpos[m->jnt_qposadr[jID]];
        }
    }
    // get the velocities
    mju_copy(x.data()+m->nv, d->qvel, m->nv);

    return x;
}

// getBounds: gets bounds on joint positions, actuator forces from the model
void CitoControl::getBounds()
{
    for( int i=0; i<m->nv; i++ )
    {
        int jID =  m->dof_jntid[i];
        if( m->jnt_limited[jID] )
        {
            isJFree[i] = 0;
            qposLB[i]  = m->jnt_range[jID*2];
            qposUB[i]  = m->jnt_range[jID*2+1];
        }
        else
        {
            isJFree[i] = 1;
            qposLB[i]  = 0;  // to be replaced by -infBnd in the initial guess
            qposUB[i]  = 0;  // to be replaced by +infBnd in the initial guess
        }
    }
    for( int i=0; i<m->nu; i++ )
    {
        if( m->actuator_ctrllimited[i] )
        {
            isAFree[i] = 0;
            tauLB[i]   = m->actuator_ctrlrange[i*2];
            tauUB[i]   = m->actuator_ctrlrange[i*2+1];
        }
        else
        {
            isAFree[i] = 1;
            tauLB[i]   = 0; // to be replaced by -infBnd in the initial guess
            tauUB[i]   = 0; // to be replaced by +infBnd in the initial guess
        }
    }
}