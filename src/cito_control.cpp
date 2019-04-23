// ***** DESCRIPTION ***********************************************************
// CitoControl class defines functions for calculating and setting control- and
// contact-related variables, i.e., joint torques and external forces on the free
// bodies.

#include "cito_control.h"


// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoControl::CitoControl(const mjModel* model) : m(model), cp(model), sl(model)
{
    // read contact model parameters
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    alpha = vscm["alpha"].as<double>();
    phiR  = vscm["phiR"].as<double>();
    // bound variables
    qposLB  = new double[m->nv];   qposUB  = new double[m->nv];
    tauLB   = new double[m->nu];   tauUB   = new double[m->nu];
    isJFree = new int[m->nv];      isAFree = new int[m->nu];
    // initialize Eigen variables
    h.resize(6*cp.nFree,1); hCon.resize(6*cp.nFree,1);
    x.resize(cp.n);
    pSR.setZero(); pSE.setZero(); nCS.setZero(); vRE.setZero();
}
CitoControl::~CitoControl()
{
    // delete bound variables
    delete []qposLB;        delete []qposUB;
    delete []tauLB;         delete []tauUB;
    delete []isJFree;       delete []isAFree;
}

// ***** FUNCTIONS *************************************************************
// takeStep: takes a full control step given a control input
void CitoControl::takeStep(mjData*d, const eigVd u, bool save, double compensateBias)
{
    if( save ) { sl.writeData(d); }
    for( int i=0; i<cp.ndpc; i++ )
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
      d->ctrl[i] = u(i) + compensateBias*d->qfrc_bias[cp.dAct[i]];
    }
    // contact model
    hCon.setZero();
    hCon = this->contactModel(d, u);
    // set external forces on the free bodies
    for( int i=0; i<cp.nFree; i++ )
    {
      for( int j=0; j<6; j++ )
      {
        d->xfrc_applied[cp.bFree[i]*6+j] = hCon(i*6+j);
      }
    }
}

// contactModel: returns contact wrench given current state and control input
eigMd CitoControl::contactModel(const mjData* d, const eigVd u)
{
    h.setZero();
    // loop for each contact pair
    for( int pI=0; pI<cp.nPair; pI++ )
    {
        // calculate the distance between the sites
        phi = getPDistance(d, pI);
        // normal force in the contact frame
        gamma = u(m->nu+pI)*exp(-alpha*phi);
        // contact generalized in the world frame
        lambda = gamma*nCS;
        // loop for each free body
        for( int fI=0; fI<cp.nFree; fI++ )
        {
          for( int i=0; i<3; i++ )
          {
            pBF[i] = d->qpos[cp.pFree[fI]+i];          // position of the center of mass of the free body
          }
          vEF = pBF - pSR;                             // vector from the end effector to the free body
          // wrench on the free body due to the contact pI: [lambda; cross(vEF, lambda)]
          h(fI*6+0) += lambda[0];
          h(fI*6+1) += lambda[1];
          h(fI*6+2) += lambda[2];
          h(fI*6+3) += -vEF[2]*lambda[1] + vEF[1]*lambda[2];
          h(fI*6+4) +=  vEF[2]*lambda[0] - vEF[0]*lambda[2];
          h(fI*6+5) += -vEF[1]*lambda[0] + vEF[0]*lambda[1];
        }
    }
    return h;
}

// getPDistance: returns the planar distance between sites given data and a pair ID
double CitoControl::getPDistance(const mjData* d, int pairID)
{
    // vectors in the world frame
    for( int i=0; i<2; i++ )
    {
        pSR[i] = d->site_xpos[cp.sPair1[pairID]*3+i];   // position of the site on the robot
        pSE[i] = d->site_xpos[cp.sPair2[pairID]*3+i];   // position of the site in the environment
        nCS[i] = cp.nCS.col(pairID)[i];                 // contact surface normal
    }
    nCS.normalize();
    // vector from the site on the robot to the site in the environment
    vRE  = pSE - pSR;
    // distance
    phiE = vRE.norm();                                  // Euclidean distance between the end effector and the environment
    phiN = vRE.dot(nCS);                                // normal distance between the end effector and the environment
    zeta  = tanh(phiR*phiE);                            // semi-sphere based on the Euclidean distance
    phiC = zeta*phiE + (1-zeta)*phiN;                   // combined distance
    // subtract the base radius of HSR
    phiC = phiC - 0.215;
    return phiC;
}

// getDistance: returns the distance between sites given data and a pair ID
double CitoControl::getDistance(const mjData* d, int pairID)
{
    // vectors in the world frame
    for( int i=0; i<3; i++ )
    {
        pSR[i] = d->site_xpos[cp.sPair1[pairID]*3+i];   // position of the site on the robot
        pSE[i] = d->site_xpos[cp.sPair2[pairID]*3+i];   // position of the site in the environment
        nCS[i] = cp.nCS.col(pairID)[i];                 // contact surface normal
    }
    // vector from the site on the robot to the site in the environment
    vRE  = pSE - pSR;
    // distance
    phiE = vRE.norm();                                  // Euclidean distance between the end effector and the environment
    phiN = vRE.dot(nCS);                                // normal distance between the end effector and the environment
    zeta  = tanh(phiR*phiE);                            // semi-sphere based on the Euclidean distance
    phiC = zeta*phiE + (1-zeta)*phiN;                   // combined distance
    return phiC;
}

// getState: converts free joints' quaternions to Euler angles so that
// the dimensionality of the state vector is 2*nv instead of nq+nv
eigVm CitoControl::getState(const mjData* d)
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
            mju_copy(x.block<3,1>(i,0).data(), d->qpos+i, 3);
            mju_copy(jFreeQuat.data(), d->qpos+i+3, 4);
            // calculate the Euler angles from the quaternion
            x(i+3) = atan2(2*(jFreeQuat[0]*jFreeQuat[1]+jFreeQuat[2]*jFreeQuat[3]), 1-2*(pow(jFreeQuat[1],2)+pow(jFreeQuat[2],2)));
            x(i+4) =  asin(2*(jFreeQuat[0]*jFreeQuat[2]-jFreeQuat[3]*jFreeQuat[1]));
            x(i+5) = atan2(2*(jFreeQuat[0]*jFreeQuat[3]+jFreeQuat[1]*jFreeQuat[2]), 1-2*(pow(jFreeQuat[2],2)+pow(jFreeQuat[3],2)));
            i += 6;             // proceed to next joint
            freeNo++;           // free joint counter
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