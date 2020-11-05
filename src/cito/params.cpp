// ***** DESCRIPTION ***********************************************************
// Params parses the model and config files and defines parameters and utility
// functions that are used across the CITO classes.

#include "cito/params.h"

// ***** CONSTRUCTOR ***********************************************************
Params::Params(const mjModel *model_) : model(model_)
{
    // read config files
    YAML::Node params = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/params.yaml");
    // task parameters
    taskType = params["taskType"].as<int>();
    desiredPos.resize(6);
    std::vector<double> desiredPosInput = {params["desiredFinalPos"].as<std::vector<double>>()};
    desiredPos = Eigen::Map<Eigen::VectorXd>(desiredPosInput.data(), desiredPosInput.size());
    controlJointDOF0 = params["controlJointDOF0"].as<int>();
    weight[0] = params["w1"].as<double>();
    weight[1] = params["w2"].as<double>();
    weight[2] = params["w3"].as<double>();
    weight[3] = params["w4"].as<double>();
    // simulation parameters
    tf = params["tf"].as<double>();
    tc = params["tc"].as<double>();
    dt = model->opt.timestep;
    N = (int)floor(tf / tc);
    ndpc = (int)floor(tc / dt);
    // read model parameters
    nu = model->nu;
    nv = model->nv;
    quatAdr = new int[nv];
    dofAdr = new int[nv];
    // get the indices of free bodies and the addresses of positions and DOF of free joints
    nFree = 0;
    for (int i = 0; i < model->nv; i++)
    {
        int jID = model->dof_jntid[i];
        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        quatAdr[i] = -1;
        dofAdr[i] = 0;
        if (model->jnt_type[jID] == mjJNT_FREE && i >= model->jnt_dofadr[jID] + 3)
        {
            if (i == model->jnt_dofadr[jID] + 3)
            {
                nFree++;
            }
            quatAdr[i] = model->jnt_qposadr[jID] + 3;
            dofAdr[i] = i - model->jnt_dofadr[jID] - 3;
        }
    }
    bFree = new int[nFree]; // indices of free bodies
    pFree = new int[nFree]; // position addresses of free DOF
    int iFree = 0;
    for (int i = 0; i < nv; i++)
    {
        int jID = model->dof_jntid[i];
        if (model->jnt_type[jID] == mjJNT_FREE)
        {
            pFree[iFree] = model->jnt_qposadr[jID];
            bFree[iFree] = model->jnt_bodyid[jID];
            iFree++;
            i += 5;
        }
    }
    // get the indices of actuated DOF
    dAct = new int[nu]; // velocity addresses of actuated DOF
    pAct = new int[nu]; // position addresses of actuated DOF
    for (int i = 0; i < nu; i++)
    {
        dAct[i] = model->jnt_dofadr[model->actuator_trnid[i * 2]];
        pAct[i] = model->jnt_qposadr[model->actuator_trnid[i * 2]];
    }
    // contact model parameters
    nPair = params["npair"].as<int>(); // number of contact pairs
    std::vector<int> site_rbt = {params["spair1"].as<std::vector<int>>()};
    std::vector<int> site_env = {params["spair2"].as<std::vector<int>>()};
    for (int i = 0; i < nPair; i++)
    {
        Eigen::Vector2i site_pair(site_rbt[i], site_env[i]);
        sites.push_back(site_pair);
    }
    // default contact surface normals in the environment
    nCS0.resize(3, nPair); // contact surface normals
    mjtNum surface_normal[3] = {1, 0, 0};
    mjtNum site_quat[4];
    for (int i = 0; i < nPair; i++)
    {
        mju_copy4(site_quat, model->site_quat + 4 * sites[i][1]);
        mju_rotVecQuat(nCS0.col(i).data(), surface_normal, site_quat);
    }
    // dimensions
    n = 2 * nv;                  // dimensionality of states
    m = nu + nPair;              // dimensionality of controls
    nTraj = (N + 1) * n + N * m; // number of parameters in discretized trajectory
}
// ***** DESTRUCTOR ************************************************************
Params::~Params()
{
    delete[] quatAdr;
    delete[] dofAdr;
    delete[] bFree;
    delete[] pFree;
    delete[] dAct;
}

// Utility functions
// skew: returns the skew symmetric matrix representation of a 3D vector
Eigen::Matrix3d Params::skew(const Eigen::Vector3d &a)
{
    Eigen::Matrix3d Ahat;
    Ahat.setZero();
    Ahat(0, 1) = -a[2];
    Ahat(0, 2) = a[1];
    Ahat(1, 2) = -a[0];
    Ahat(1, 0) = a[2];
    Ahat(2, 0) = -a[1];
    Ahat(2, 1) = a[0];
    return Ahat;
}

// skewCross: performs and returns a x b using skew-symmetric transformation
Eigen::Vector3d Params::skewCross(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
{
    Eigen::Vector3d c;
    c[0] = -a[2] * b[1] + a[1] * b[2];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = -a[1] * b[0] + a[0] * b[1];
    return c;
}

// quat2Euler: converts a quaternion (w,x,y,z) into ZYX Euler angles
Eigen::Vector3d Params::quat2Euler(const Eigen::Vector4d &q)
{
    Eigen::Vector3d e;
    e[0] = atan2(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (pow(q[1], 2) + pow(q[2], 2))); // roll
    e[1] = asin(2 * (q[0] * q[2] - q[3] * q[1]));                                         // pitch
    e[2] = atan2(2 * (q[0] * q[3] + q[1] * q[2]), 1 - 2 * (pow(q[2], 2) + pow(q[3], 2))); // yaw
    return e;
}

// evalNormalJac: calculates the contact normal Jacobian w.r.t. rotational DOF
Eigen::Matrix3d Params::evalNormalJac(const Eigen::Vector4d &q, int pair)
{
    Eigen::Vector3d e = quat2Euler(q); // roll, pitch, yaw
    Eigen::Matrix3d Rx, Ry, Rz, dRx, dR_dx, dR_dw;
    dRx.setZero();
    dRx(1, 1) = -sin(e[0]);
    dRx(1, 2) = -cos(e[0]);
    dRx(2, 1) = cos(e[0]);
    dRx(2, 2) = -sin(e[0]);
    Ry.setZero();
    Ry(0, 0) = cos(e[1]);
    Ry(0, 2) = sin(e[1]);
    Ry(1, 1) = 1.;
    Ry(2, 0) = -sin(e[1]);
    Ry(2, 2) = cos(e[1]);
    Rz.setZero();
    Rz(0, 0) = cos(e[2]);
    Rz(0, 1) = -sin(e[2]);
    Rz(1, 0) = sin(e[2]);
    Rz(1, 1) = cos(e[2]);
    Rz(2, 2) = 1.;
    // calculate only the x component of the tensor since it's multiplied by ux
    dR_dx = Rz * Ry * dRx;
    // not sure why there is the negative sign but this matches the num. diff.
    dR_dw.setZero();
    dR_dw.col(2) = -dR_dx.col(2);
    // This assumes site orientations align with a world axis in the model
    // and is tested for only +/- x/y directions and not for z axis.
    if (abs(nCS0(0, pair)) < 1e-6)
        dR_dw.col(0) = nCS0(1, pair) * dR_dx.col(1);
    else
        dR_dw.col(1) = -nCS0(0, pair) * dR_dx.col(1);
    return dR_dw;
}

// getState: converts free joints' quaternions to Euler angles so that
// the dimensionality of the state vector is 2*nv instead of nq+nv
eigVd Params::getState(const mjData *d)
{
    // create the state and quaternion vectors
    eigVd x(n), jFreeQuat(4);
    x.setZero();
    // get the positions
    int freeNo = 0;
    for (int i = 0; i < model->nq; i++)
    {
        int jID = model->dof_jntid[i - freeNo];
        if (model->jnt_type[jID] == mjJNT_FREE)
        {
            mju_copy3(x.segment(i, 3).data(), d->qpos + i);
            mju_copy4(jFreeQuat.data(), d->qpos + i + 3);
            // convert quaternion into Euler angles
            x.segment(i + 3, 3) = quat2Euler(jFreeQuat);
            // count the free joint and proceed to next joint
            freeNo++;
            i += 6;
        }
        else
        {
            x(i - freeNo) = d->qpos[model->jnt_qposadr[jID]];
        }
    }
    // get the velocities
    mju_copy(x.data() + model->nv, d->qvel, model->nv);
    // return the state w/ Euler angles
    return x;
}