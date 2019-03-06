// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include "cito_scvx.h"

// ***** MuJoCo model & data ****************************************************/
mjModel *m = NULL;
mjData  *d = NULL;
//================================================================================
int main(int argc, char const *argv[]) {
    // ***** Model file **********************************************************/
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/"  + paths::modelFile;
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";
    // ***** Trajectories ********************************************************/
    ctrlTraj U0, U; U0.resize(NTS); U.resize(NTS);
    trajectory traj;
    // ***** MuJoCo initialization ***********************************************/
    // Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Load xml model
    if( strlen(modelPath)>4 && !strcmp(modelPath+strlen(modelPath)-4, ".mjb") )
    {       m = mj_loadModel(modelPath, NULL); }
    else {  m = mj_loadXML(modelPath, NULL, NULL, 0); }
    if( !m ) { mju_error("Cannot load the model"); }
    // Create data
    d = mj_makeData(m);
    // ***** Create objects for CITO *********************************************/
    CitoSCvx scvx(m);
    // ***** Initial control trajectory ******************************************/
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    double kCon0 = vscm["kCon0"].as<double>();
    for (int i = 0; i < NTS; i++)
    {
        U0[i].setZero();
        for (int j = 0; j < NPAIR; j++)
        {
            U0[i][NU + j] = kCon0;
        }
    }
    // ***** Run successive convexification **************************************/
    U = scvx.solveSCvx(U0);
    // ***** Evaluate the optimal trajectory *************************************/
    traj = scvx.runSimulation(U, false, true);
    // Print the trajectory
    std::cout << "\n\nOptimal trajectory:\n";
    for( int i=0; i<NTS; i++ )
    {
        std::cout << "time step " << i << ":\n\t\tpos = " << traj.X[i].block<NV, 1>(0, 0).transpose() << "\n";
        std::cout << "\t\t vel = " << traj.X[i].block<NV, 1>(NV, 0).transpose() << "\n";
        std::cout << "\t\t tau = ";
        std::cout << traj.U[i].block<NU, 1>(0, 0).transpose() << "\n";
        std::cout <<"\t\t KCon = ";
        std::cout << traj.U[i].block<NPAIR, 1>(NU, 0).transpose() << "\n\n";
    }
    std::cout << "time step " << NTS << ":\n\t\tpos = " << traj.X[NTS].block<NV, 1>(0, 0).transpose() << "\n";
    std::cout << "\t\t vel = " << traj.X[NTS].block<NV, 1>(NV, 0).transpose() << "\n";
    // ***** Evaluate the optimal const ******************************************/
    double J = scvx.getCost(traj.X[NTS], traj.U);
    std::cout << "J = " << J;
    // ***** MuJoCo shut down ****************************************************/
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
