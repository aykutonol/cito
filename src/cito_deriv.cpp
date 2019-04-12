#include "cito_deriv.h"

// ***** DESCRIPTION ***********************************************************
// CitoDeriv class defines functions for calculating the derivatives of the
// dynamics using Pinocchio.

// ***** CONSTRUCTOR ***********************************************************
CitoDeriv::CitoDeriv(const mjModel* mModel) : m(mModel), cp(mModel)
{
    // ***** Pinocchio initialization ********************************************/
    YAML::Node params = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/params.yaml");
    std::string urdfPathStr = paths::workspaceDir + "/src/cito/model/" + params["urdf"].as<std::string>();
    const char *urdfPath = urdfPathStr.c_str();
    std::cout << "\n\nURDF path for Pinocchio: " << urdfPath << "\n\n\n";
    pinocchio::urdf::buildModel(urdfPath, model);
    // time coefficient
    tM = cp.dt*cp.ndpc;
    // resize Eigen variables
    q.resize(m->nq); v.resize(m->nv); tau.resize(m->nu);
    Fx.resize(cp.n, cp.n); Fu.resize(cp.n, cp.m);
}

// ***** FUNCTIONS *************************************************************
// linDyn: calculates derivatives of the state and control trajectories
void CitoDeriv::linDyn(const mjData* d, const eigVd u, mjtNum* Fxd, mjtNum* Fud, double compensateBias)
{
    pinocchio::Data data(model);
    mju_copy(q.data(), d->qpos, m->nq);
    mju_copy(v.data(), d->qvel, m->nv);
    for( int i=0; i<m->nu; i++ )
    {
        tau(i) = u(i) + compensateBias*d->qfrc_bias[cp.dAct[i]];
    }
    pinocchio::computeABADerivatives(model, data, q, v, tau);
    // initialize derivative matrices
    Fx.setZero();   Fu.setZero();
    // build state derivative matrix
    Fx.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
    Fx.topRightCorner(m->nv, m->nv)    = tM*Eigen::MatrixXd::Identity(m->nv, m->nv);
    Fx.bottomLeftCorner(m->nv, m->nv)  = tM*data.ddq_dq;
    Fx.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + tM*data.ddq_dv;
    // build control derivative matrix
    Fu.bottomRows(m->nv) = tM*data.Minv;
    // return matrices
    mju_copy(Fxd, Fx.data(), cp.n*cp.n);
    mju_copy(Fud, Fu.data(), cp.n*cp.m);
}