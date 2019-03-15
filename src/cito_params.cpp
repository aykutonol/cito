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
    // contact model parameters
    npair = vscm["npair"].as<int>();
    // dimensions
    n = 2*NV;               // dimensionality of states
    m = nu+npair;           // dimensionality of controls
    ntraj = (N+1)*n + N*m;  // number of parameters in discretized trajectory
}