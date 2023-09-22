/*! Main */
/**
 *  \brief main solves the CITO problem using SCVX
 *
 *  This file initializes MuJoCo and runs the successive convexification algorithm.
 *  The optimal trajectory is recorded to cito_ws/logs for playback and execution.
 *
 *  \author Aykut Onol
 */

#include "cito/penalty_loop.h"

int main(int argc, char const *argv[])
{
    int testNumber = -1;
    std::string keypoint_method = "SI2";
    for(int i = 1; i < argc; i++){
        testNumber = std::stoi(argv[i]);
    }
    std::cout << "test number: " << testNumber << "\n";

    // ***** MuJoCo initialization ***********************************************/
    // Model file
    YAML::Node params = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/params.yaml");
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/" + params["model"].as<std::string>();
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";
    // Load the model
    mjModel *m = NULL;
    if (strlen(modelPath) > 4 && !strcmp(modelPath + strlen(modelPath) - 4, ".mjb"))
    {
        m = mj_loadModel(modelPath, NULL);
    }
    else
    {
        m = mj_loadXML(modelPath, NULL, NULL, 0);
    }
    if (!m)
    {
        mju_error("Cannot load the model");
    }
    // ***** Create objects for CITO *********************************************/
    Params cp(m, testNumber);
    Control cc(m, &cp);
    SCVX scvx(m, &cp, &cc);
    PenaltyLoop pl(m, &cp, &cc, &scvx);
    // ***** Trajectories ********************************************************/
    eigMd U0, U;
    U0.resize(cp.m, cp.N);
    U.resize(cp.m, cp.N);
    trajectory traj;
    // ***** Initial control trajectory ******************************************/
    YAML::Node vscm = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/vscm.yaml");
    double kCon0 = vscm["kCon0"].as<double>();
    for (int i = 0; i < cp.N; i++)
    {
        U0.col(i).setZero();
        for (int j = 0; j < cp.nPair; j++)
        {
            U0.col(i)[m->nu + j] = kCon0;
        }
    }
    // ***** Run the planner *****************************************************/
    auto tPlanStart = std::chrono::system_clock::now();
    auto cpuPlanStart = std::clock();
    if (params["penaltyLoop"].as<bool>())
        U = pl.solve(U0);
    else
        U = scvx.solveSCVX(U0);
    auto cpuPlanEnd = std::clock();
    auto tPlanEnd = std::chrono::system_clock::now();
    // ***** Evaluate the optimal trajectory *************************************/
    traj = scvx.runSimulation(U, false, 2, 1);
    // Print the trajectory
//    std::cout << "\n\nOptimal trajectory:\n";
//    for (int i = 0; i < cp.N; i++)
//    {
//        std::cout << "time step " << i << ":\n\t\tpos = " << traj.X.col(i).head(m->nv).transpose() << "\n";
//        std::cout << "\t\t vel = " << traj.X.col(i).tail(m->nv).transpose() << "\n";
//        std::cout << "\t\t tau = ";
//        std::cout << traj.U.col(i).head(m->nu).transpose() << "\n";
//        std::cout << "\t\t KCon = ";
//        std::cout << traj.U.col(i).tail(cp.nPair).transpose() << "\n\n";
//    }
//    std::cout << "time step " << cp.N << ":\n\t\tpos = " << traj.X.col(cp.N).head(m->nv).transpose() << "\n";
//    std::cout << "\t\t vel = " << traj.X.col(cp.N).tail(m->nv).transpose() << "\n";
    // ***** Evaluate the optimal const ******************************************/
    double J = scvx.getCost(traj.X, traj.U);
    std::cout << "\nJ = " << J << ", kmax = " << U.bottomRows(cp.nPair).maxCoeff() << "\n\nINFO: Planning completed in " << std::chrono::duration<double>(tPlanEnd - tPlanStart).count() << " wall-clock s and " << 1000.0 * (cpuPlanEnd - cpuPlanStart) / CLOCKS_PER_SEC << " CPU ms.\n\n";
    // ***** MuJoCo shut down ****************************************************/
    mj_deleteModel(m);

    // Save the testing data to the correct folder
    double costReduction = scvx.costReduction;
    double optTime = scvx.optTime;
    double derivsTime = scvx.derivsTime;
    double qpTime = scvx.qpTime;

    if(testNumber != -1){
        std::string testDir = paths::workspaceDir + "/src/cito/testingData/data/" + keypoint_method + "/";

        // Open a csv file and save the data
        std::ofstream csvFile;
        std::string file_name = testDir + std::to_string(testNumber) + ".csv";
        std::cout << "filename: " << file_name << "\n";
        csvFile.open(testDir + std::to_string(testNumber) + ".csv");

        // Save the data - cost reduction, opt time, derivs time, qp time
        csvFile << optTime << ",";
        csvFile << costReduction << ",";
        csvFile << derivsTime << ",";
        csvFile << qpTime << "\n";

    }


    return 0;
}
