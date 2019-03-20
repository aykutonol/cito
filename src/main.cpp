/*! Main */
/**
 *  \brief main solves the CITO problem
 *
 *  This file initializes MuJoCo and runs the successive convexification algorithm.
 *  The optimal trajectory is recorded to cito_ws/logs for playback and execution.
 *
 *  \author Aykut Onol
 */

#include "cito_scvx.h"

// MuJoCo model
mjModel *m = NULL;

int main(int argc, char const *argv[]) {
    // ***** MuJoCo initialization ***********************************************/
    // Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Model file
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/"  + paths::modelFile;
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";
    // Load the model
    if( strlen(modelPath)>4 && !strcmp(modelPath+strlen(modelPath)-4, ".mjb") )
    {       m = mj_loadModel(modelPath, NULL); }
    else {  m = mj_loadXML(modelPath, NULL, NULL, 0); }
    if( !m ) { mju_error("Cannot load the model"); }
    // ***** Create objects for CITO *********************************************/
    CitoParams cp(m);
    CitoSCvx scvx(m);
    // ***** Trajectories ********************************************************/
    eigDbl U0, U; U0.resize(cp.m,cp.N); U.resize(cp.m,cp.N);
    trajectory traj;
    // ***** Initial control trajectory ******************************************/
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/vscm.yaml");
    double kCon0 = vscm["kCon0"].as<double>();
    for (int i=0; i<cp.N; i++)
    {
        U0.col(i).setZero();
        for (int j=0; j<cp.nPair; j++)
        {
            U0.col(i)[m->nu + j] = kCon0;
        }
    }
    // ***** Run successive convexification **************************************/
    U = scvx.solveSCvx(U0);
    // ***** Evaluate the optimal trajectory *************************************/
    traj = scvx.runSimulation(U, false, true);
    // Print the trajectory
    std::cout << "\n\nOptimal trajectory:\n";
    for( int i=0; i<cp.N; i++ )
    {
        std::cout << "time step " << i << ":\n\t\tpos = " << traj.X.col(i).block<NV,1>(0,0).transpose() << "\n";
        std::cout << "\t\t vel = " << traj.X.col(i).block<NV, 1>(m->nv,0).transpose() << "\n";
        std::cout << "\t\t tau = ";
        std::cout << traj.U.col(i).block<NU,1>(0,0).transpose() << "\n";
        std::cout <<"\t\t KCon = ";
        std::cout << traj.U.col(i).block<NPAIR,1>(m->nu,0).transpose() << "\n\n";
    }
    std::cout << "time step " << cp.N << ":\n\t\tpos = " << traj.X.col(cp.N).block<NV,1>(0,0).transpose() << "\n";
    std::cout << "\t\t vel = " << traj.X.col(cp.N).block<NV, 1>(m->nv, 0).transpose() << "\n";
    // ***** Evaluate the optimal const ******************************************/
    double J = scvx.getCost(traj.X.col(cp.N), traj.U);
    std::cout << "J = " << J << "\n\nINFO: Planning completed.\n\n";
    // ***** MuJoCo shut down ****************************************************/
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
