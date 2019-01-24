// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include <iostream>

#include "cito_scvx.h"

// ********* mujoco model & data ************************************************/
mjModel *m = NULL;
mjData  *d = NULL;
// ********* state, control, and derivative threads *****************************/
stateVecThread X, XTilde;       ctrlVecThread U0, Uopt;
stateDerThread Fx;          ctrlDerThread Fu;
//==============================================================================
int main(int argc, char const *argv[]) {
    // ********* mujoco initialization ******************************************/
    // activate mujoco
    mj_activate(paths::mjKey);
    // load xml model
    if( strlen(paths::modelFile)>4 && !strcmp(paths::modelFile+strlen(paths::modelFile)-4, ".mjb") )
    {       m = mj_loadModel(paths::modelFile, NULL); }
    else {  m = mj_loadXML(paths::modelFile, NULL, NULL, 0); }
    if( !m ) { mju_error("Cannot load the model"); }
    // create data
    d = mj_makeData(m);
    // ********* create objects for CITO ****************************************/
    CitoControl cc(m);
    CitoNumDiff nd(m);
    CitoSCvx    scvx(m);
    // create trajectories
    U0.resize(NTS); Uopt.resize(NTS); X.resize(NTS+1); XTilde.resize(NTS+1);
    Fx.resize(NTS); Fu.resize(NTS);
    trajectory trajOpt;
    // ********* initial control trajectory *************************************/
    for (int i = 0; i < NTS; i++)
    {
        U0[i].setZero();
        for (int j = 0; j < NPAIR; j++)
        {
            U0[i][NU + j] = params::kCon0[j];
        }
    }
    // ********* simulation *****************************************************/
//    for (int i = 0; i < NTS; i++)
//    {
//        for (int j = 0; j < params::ndpc; j++)
//        {
//            mj_step1(m, d);
//            cc.setControl(d, U0[i]);
//            mj_step2(m, d);
//        }
//        X[i].setZero(); XTilde[i].setZero(); Fx[i].setZero(); Fu[i].setZero();
//        X[i] = cc.getState(d);
//        std::cout << "i: " << i << "\t tau: ";
//        mju_printMat(d->ctrl, 1, m->nu);
//        std::cout <<"\t\t kCon: ";
//        std::cout << U0[i].block<NPAIR, 1>(NU, 0).transpose() << "\n\n";
//        std::cout << "\t\t xfrc: ";
//        mju_printMat(d->xfrc_applied + params::bfree[0] * 6, 1, 6);
//        std::cout << "\t\t torso pose: ";
//        mju_printMat(d->qpos + params::jfree[0], 1, 7);
//        std::cout << "X: " << X[i].transpose() << "\n\n";
////        nd.linDyn(d, U[i], Fx[i].data(), Fu[i].data());
////        std::cout << "Fx:\n" << Fx[i] << '\n';
////        std::cout << "Fu:\n" << Fu[i] << '\n';
//    }
    Uopt = scvx.solveSCvx(U0);
    trajOpt = scvx.runSimulation(Uopt, false, true);
    std::cout << "\n\nOptimal trajectory:\n";
    for( int i=0; i<NTS; i++ )
    {
        std::cout << "time step" << i << ":\n\tX = " << trajOpt.X[i].transpose() << "\n";
        std::cout << "\t\t tau = ";
        std::cout << trajOpt.U[i].block<NU, 1>(0, 0).transpose() << "\n";
        std::cout <<"\t\t kCon = ";
        std::cout << trajOpt.U[i].block<NPAIR, 1>(NU, 0).transpose() << "\n\n";
    }
    std::cout << "X" << NTS << ": " << trajOpt.X[NTS].transpose() << "\n\n";

    double J = scvx.getCost(trajOpt.X[NTS], trajOpt.U);
    std::cout << "J = " << J;
    std::cout << "\n\nNTS: " << NTS << ", N: " << N << ", M: " << M << ", NU: " << NU << ", NPAIR: " << NPAIR << "\n";
    // ********* mujoco shut down ***********************************************/
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
