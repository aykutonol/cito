#include <gtest/gtest.h>

#include "cito_scvx.h"

// MuJoCo model and data
mjModel *m = NULL;
mjData  *d = NULL;
int mj_activated = 0;

// YAML file for parameters
YAML::Node params;

// Assert model and data are not empty
TEST(CitoTests, check_mujoco) {
    ASSERT_TRUE(mj_activated>0);
    ASSERT_TRUE(m!=NULL);
    ASSERT_TRUE(d!=NULL);
}

// Assert params file loaded
TEST(CitoTests, check_params) {
    ASSERT_FALSE(params.IsNull());
}

int main(int argc, char **argv) {
    // Initialize testing
    testing::InitGoogleTest(&argc, argv);

    // Model file
    params = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/params.yaml");
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/" + params["model"].as<std::string>();
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";

    // Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activated = mj_activate(mjKeyPath);
    // Load the model
    if( strlen(modelPath)>4 && !strcmp(modelPath+strlen(modelPath)-4, ".mjb") )
    {       m = mj_loadModel(modelPath, NULL); }
    else {  m = mj_loadXML(modelPath, NULL, NULL, 0); }
    if( !m ) { mju_error("Cannot load the model"); }
    // Create data
    d = mj_makeData(m);

    // Run and return tests
    return RUN_ALL_TESTS();
}