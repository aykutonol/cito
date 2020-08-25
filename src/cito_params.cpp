// ***** DESCRIPTION ***********************************************************
// CitoParams class parses the model and config files and defines parameters
// that are used across classes.

#include "cito_params.h"

// ***** CONSTRUCTOR ***********************************************************
CitoParams::CitoParams(const mjModel* model_) : model(model_)
{
    // read config files
    YAML::Node params = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/params.yaml");
    // task parameters
    desiredPos.resize(6);
    std::vector<double> desiredPosInput = { params["desiredFinalPos"].as<std::vector<double>>() };
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
    N    = (int) floor(tf/tc);
    ndpc = (int) floor(tc/dt);
    // read model parameters
    nu = model->nu;
    nv = model->nv;
    quatAdr = new int[nv];
    dofAdr  = new int[nv];
    // get the indices of free bodies and the addresses of positions and DOF of free joints
    nFree = 0;
    for( int i=0; i<model->nv; i++ )
    {
        int jID = model->dof_jntid[i];
        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        quatAdr[i] = -1; dofAdr[i] = 0;
        if( model->jnt_type[jID]==mjJNT_FREE && i>=model->jnt_dofadr[jID]+3 )
        {
            if( i== model->jnt_dofadr[jID]+3 )
            {
                nFree++;
            }
            quatAdr[i] = model->jnt_qposadr[jID] + 3;
            dofAdr[i]  = i - model->jnt_dofadr[jID] - 3;
        }
    }
    bFree = new int[nFree];     // indices of free bodies
    pFree = new int[nFree];     // position addresses of free DOF
    int iFree = 0;
    for( int i=0; i<nv; i++ )
    {
        int jID = model->dof_jntid[i];
        if( model->jnt_type[jID]==mjJNT_FREE )
        {
            pFree[iFree] = model->jnt_qposadr[jID];
            bFree[iFree] = model->jnt_bodyid[jID];
            iFree++;
            i += 5;
        }
    }
    // get the indices of actuated DOF
    dAct = new int[nu];         // addresses of actuated DOF
    for( int i=0; i<nu; i++ )
    {
        dAct[i] = model->jnt_dofadr[model->actuator_trnid[i*2]];
    }
    // contact model parameters
    nPair  = params["npair"].as<int>();     // number of contact pairs
    std::vector<int> site_rbt = { params["spair1"].as<std::vector<int>>() };
    std::vector<int> site_env = { params["spair2"].as<std::vector<int>>() };
    for(int i=0; i<nPair; i++)
    {
        Eigen::Vector2i site_pair(site_rbt[i], site_env[i]);
        sites.push_back(site_pair);
    }
    // contact surface normals in the environment
    nCS.resize(3, nPair);                   // contact surface normals
    mjtNum surface_normal[3] = {1, 0, 0};
    mjtNum site_quat[4];
    for( int i=0; i<nPair; i++ )
    {
        mju_copy4(site_quat, model->site_quat+4*sites[i][1]);
        mju_rotVecQuat(nCS.col(i).data(), surface_normal, site_quat);
    }
    // dimensions
    n = 2*nv;               // dimensionality of states
    m = nu+nPair;           // dimensionality of controls
    nTraj = (N+1)*n + N*m;  // number of parameters in discretized trajectory
}
// ***** DESTRUCTOR ************************************************************
CitoParams::~CitoParams()
{
    delete[] quatAdr;   delete[] dofAdr;
    delete[] bFree;     delete[] pFree;
    delete[] dAct;
}

// Utility functions
// skewCross: performs and returns a x b using skew-symmetric transformation
Eigen::Vector3d CitoParams::skewCross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    Eigen::Vector3d c;
    c[0] = -a[2]*b[1] + a[1]*b[2];
    c[1] =  a[2]*b[0] - a[0]*b[2];
    c[2] = -a[1]*b[0] + a[0]*b[1];
    return c;
}

// quat2Euler: converts a quaternion (w,x,y,z) into Euler angles
Eigen::Vector3d CitoParams::quat2Euler(Eigen::Vector4d q) {
    Eigen::Vector3d e;
    e[0] = atan2(2*(q[0]*q[1]+q[2]*q[3]), 1-2*(pow(q[1],2)+pow(q[2],2)));
    e[1] =  asin(2*(q[0]*q[2]-q[3]*q[1]));
    e[2] = atan2(2*(q[0]*q[3]+q[1]*q[2]), 1-2*(pow(q[2],2)+pow(q[3],2)));
    return e;
}