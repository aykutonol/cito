#include "cito_params.h"

// ***** DESCRIPTION ***********************************************************
// CitoParams class parses the model and config files and defines parameters
// that are used across classes.

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoParams::CitoParams(const mjModel* model) : model(model)
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
    sPair1.resize(nPair);                   // indices of sites on the robot
    sPair2.resize(nPair);                   // incides of corresponding sites in the environment
    nCS.resize(3,nPair);                    // contact surface normals
    std::vector<int> stdVecInt = { params["spair1"].as<std::vector<int>>() };
    sPair1 = Eigen::Map<Eigen::VectorXi>(stdVecInt.data(), stdVecInt.size());
    stdVecInt = { params["spair2"].as<std::vector<int>>() };
    sPair2 = Eigen::Map<Eigen::VectorXi>(stdVecInt.data(), stdVecInt.size());
    mjtNum normal[3] = {1, 0, 0};
    mjtNum mjQuat[4];
    for( int i=0; i<nPair; i++ )
    {
        for( int j=0; j<4; j++ )
        {
            mjQuat[j] = model->site_quat[sPair2(i)*4+j];
        }
        mju_rotVecQuat(nCS.col(i).data(), normal, mjQuat);
    }
    // dimensions
    n = 2*nv;               // dimensionality of states
    m = nu+nPair;           // dimensionality of controls
    nTraj = (N+1)*n + N*m;  // number of parameters in discretized trajectory
}