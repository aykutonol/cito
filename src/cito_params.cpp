#include "cito_params.h"

// ***** DESCRIPTION ***********************************************************
// CitoParams class defines global variables that are specific to simulation,
// robot, and environment as well as general types and structures.

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoParams::CitoParams(const mjModel* model) : model(model)
{
    // read config files
    YAML::Node sim  = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/sim.yaml");
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    // simulation parameters
    tf = sim["tf"].as<double>();
    tc = sim["tc"].as<double>();
    dt = model->opt.timestep;
    N    = (int) floor(tf/tc);
    ndpc = (int) floor(tc/dt);
    // read model parameters
    nu = model->nu;
    nv = model->nv;
    quatAdr = new int[nv];
    dofAdr  = new int[nv];
    for( int i=0; i<nv; i++ )
    {
        int jID = model->dof_jntid[i];
        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        quatAdr[i] = -1; dofAdr[i] = 0;
        if( model->jnt_type[jID]==mjJNT_FREE && i>=model->jnt_dofadr[jID]+3 )
        {
            quatAdr[i] = model->jnt_qposadr[jID] + 3;
            dofAdr[i]  = i - model->jnt_dofadr[jID] - 3;
        }
    }
    // contact model parameters
    npair = vscm["npair"].as<int>();
    // dimensions
    n = 2*nv;               // dimensionality of states
    m = nu+npair;           // dimensionality of controls
    ntraj = (N+1)*n + N*m;  // number of parameters in discretized trajectory
}