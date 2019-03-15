#include "cito_params.h"

// ***** DESCRIPTION ***********************************************************
// CitoParams class defines global variables that are specific to simulation,
// robot, and environment as well as general types and structures.

// ***** CONSTRUCTOR & DESTRUCTOR **********************************************
CitoParams::CitoParams()
{
    // read parameters
    YAML::Node params = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/params.yaml");

    test = params["test"].as<int>();
}