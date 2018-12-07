// =============================== //
// *** Developed by Aykut Onol *** //
// =============================== //

#include <iostream>

#include "cito_control.h"

// ********* mujoco model & data **********************************************/
mjModel *m = NULL;
mjData *d = NULL;
//  ********* create objects for CITO *****************************************/
CitoControl cc;

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
    // ********* initial control trajectory *************************************/
    ctrlVecThread U;
    U.resize(NTS);
    ctrlVec_t u0;
    u0.setZero();
    for (int i = 0; i < NTS; i++) {
        U[i] = u0;
        for (int j = 0; j < params::npair; j++) {
            U[i][params::nact + j] = params::kcon0[j];
        }
        std::cout << "i: " << i << ", U = " << U[i].transpose() << '\n';
    }
    // ********* simulation *****************************************************/
    for (int i = 0; i < NTS; i++) {
        for (int j = 0; j < params::ndpc; j++) {
            mj_step1(m, d);
            cc.setControl(d, U[i]);
            mj_step2(m, d);
        }
        std::cout << "i: " << i << "\t tau: ";
        mju_printMat(d->ctrl, 1, m->nu);
        std::cout << "\t xfrc: ";
        mju_printMat(d->xfrc_applied + params::bfree[0] * 6, 1, 6);
        std::cout << "\t torso pose: ";
        mju_printMat(d->qpos + params::jfree[0], 1, 7);
    }
    // ********* mujoco shut down ***********************************************/
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();

    return 0;
}
