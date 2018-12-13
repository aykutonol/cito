// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include <iostream>

#include "cito_numdiff.h"

// ********* mujoco model & data ************************************************/
mjModel *m = NULL;
mjData  *d = NULL;
// ********* state, control, and derivative threads *****************************/
stateVecThread X, XL, dX;   ctrlVecThread U, dU, Utemp;
stateMatThread Fx;          ctrlMatThread Fu;
//==============================================================================
int main(int argc, char const *argv[]) {
    // ********* mujoco initialization ******************************************/
    // activate mujoco
    mj_activate("/home/aykut/Development/cito/src/cito/bin/mjkey.txt");
    // load xml model
    m = mj_loadXML("/home/aykut/Development/cito/src/cito/model/flymanoid.xml", 0, 0, 0);
    if (!m)
        mju_error("Cannot load the model");
    // create data
    d = mj_makeData(m);
    // ********* create objects for CITO ****************************************/
    CitoControl cc(m);
    CitoNumDiff nd(m);
    // create trajectories
    U.resize(NTS); X.resize(NTS+1); XL.resize(NTS+1); Fx.resize(NTS); Fu.resize(NTS);
    // ********* initial control trajectory *************************************/
    for (int i = 0; i < NTS; i++)
    {
        U[i].setZero();
        for (int j = 0; j < params::npair; j++) {
            U[i][params::nact + j] = params::kcon0[j];
        }
    }
    // ********* simulation *****************************************************/
    for (int i = 0; i < NTS; i++) {
        for (int j = 0; j < params::ndpc; j++)
        {
            mj_step1(m, d);
            cc.setControl(d, U[i]);
            mj_step2(m, d);
        }
        X[i].setZero(); XL[i].setZero(); Fx[i].setZero(); Fu[i].setZero();
        X[i] = cc.getState(d);
        std::cout << "i: " << i << "\t tau: ";
        mju_printMat(d->ctrl, 1, m->nu);
        std::cout << "\t\t xfrc: ";
        mju_printMat(d->xfrc_applied + params::bfree[0] * 6, 1, 6);
        std::cout << "\t\t torso pose: ";
        mju_printMat(d->qpos + params::jfree[0], 1, 7);
        std::cout << "X: " << X[i].transpose() << '\n';
        nd.linDyn(d, U[i], Fx[i].data(), Fu[i].data());
        std::cout << "Fx:\n" << Fx[i] << '\n';
        std::cout << "Fu:\n" << Fu[i] << '\n';
    }
    // ********* mujoco shut down ***********************************************/
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
