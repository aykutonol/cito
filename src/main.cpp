// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include "cito_scvx.h"

// ********* mujoco model & data ************************************************/
mjModel *m = NULL;
mjData  *d = NULL;
//==============================================================================
int main(int argc, char const *argv[]) {
    // ********* model file *****************************************************/
    YAML::Node paths = YAML::LoadFile(workspaceDir+"/src/cito/config/path.yaml");
    std::string modelFileTemp = paths["modelFile"].as<std::string>();
    const char *modelFile = modelFileTemp.c_str();
    std::cout << "Model file: " << modelFile << "\n";
    // ********* threads & structs **********************************************/
    ctrlVecThread U0, UOpt;
    trajectory trajOpt;
    // ********* mujoco initialization ******************************************/
    // activate mujoco
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // load xml model
    if( strlen(modelFile)>4 && !strcmp(modelFile+strlen(modelFile)-4, ".mjb") )
    {       m = mj_loadModel(modelFile, NULL); }
    else {  m = mj_loadXML(modelFile, NULL, NULL, 0); }
    if( !m ) { mju_error("Cannot load the model"); }
    // create data
    d = mj_makeData(m);
    // ********* create objects for CITO ****************************************/
//    CitoControl cc(m);
//    CitoNumDiff nd(m);
    CitoSCvx    scvx(m);
    // ********* create objects for CITO ****************************************/
    U0.resize(NTS); UOpt.resize(NTS);
    // ********* initial control trajectory *************************************/
    for (int i = 0; i < NTS; i++)
    {
        U0[i].setZero();
        for (int j = 0; j < NPAIR; j++)
        {
            U0[i][NU + j] = params::kCon0[j];
        }
    }
    // ********* run successive convexification *********************************/
    UOpt = scvx.solveSCvx(U0);
    // ********* evaluate the optimal trajectory ********************************/
    trajOpt = scvx.runSimulation(UOpt, false, true);
    // print the trajectory
    std::cout << "\n\nOptimal trajectory:\n";
    for( int i=0; i<NTS; i++ )
    {
        std::cout << "time step " << i << ":\n\t\tpos = " << trajOpt.X[i].block<NV, 1>(0, 0).transpose() << "\n";
        std::cout << "\t\t vel = " << trajOpt.X[i].block<NV, 1>(NV, 0).transpose() << "\n";
        std::cout << "\t\t tau = ";
        std::cout << trajOpt.U[i].block<NU, 1>(0, 0).transpose() << "\n";
        std::cout <<"\t\t KCon = ";
        std::cout << trajOpt.U[i].block<NPAIR, 1>(NU, 0).transpose() << "\n\n";
    }
    std::cout << "time step " << NTS << ":\n\t\tpos = " << trajOpt.X[NTS].block<NV, 1>(0, 0).transpose() << "\n";
    std::cout << "\t\t vel = " << trajOpt.X[NTS].block<NV, 1>(NV, 0).transpose() << "\n";
    // ********* evaluate the optimal const *************************************/
    double J = scvx.getCost(trajOpt.X[NTS], trajOpt.U);
    std::cout << "J = " << J;
    // ********* mujoco shut down ***********************************************/
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
