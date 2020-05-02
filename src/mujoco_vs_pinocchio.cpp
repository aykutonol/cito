#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

void printConfig(mjModel* m, mjData* d)
{
    std::cout << "\ntime: " << d->time << "\n";
    std::cout << "qpos: ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel: ";
    mju_printMat(d->qvel, 1, m->nv);
    std::cout << "ctrl: ";
    mju_printMat(d->ctrl, 1, m->nu);
    std::cout << "qacc: ";
    mju_printMat(d->qacc, 1, m->nv);
    std::cout << "qcst: ";
    mju_printMat(d->qfrc_constraint, 1, m->nv);
}

int main(int argc, char const *argv[]) {
    // MuJoCo
    /// Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    /// Load model
    std::string mjModelPathStr = paths::workspaceDir + "/src/cito/model/sawyer_push.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    /// Create data
    mjData* d = mj_makeData(m);
    /// Set the joint positions to the preset configuration in the model
        mju_copy(d->qpos, m->key_qpos, m->nq);
    /// Evaluate forward dynamics
    mj_forward(m, d);

    // Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);
    /// Create data
    Eigen::VectorXd q(model.nq), v(model.nv);
    mju_copy(q.data(), d->qpos+7, model.nq);
    v.setZero();

    // Evaluate forward dynamics and print mass matrix
    pinocchio::computeAllTerms(model, data, q, v);
    std::cout << "Pinocchio:\n" << data.M << "\n";

    // Get mass matrix from MuJoCo
    mjtNum* denseM = mj_stackAlloc(d, m->nv*m->nv);
    mj_fullM(m, denseM, d->qM);
    Eigen::MatrixXd fullM(m->nv, m->nv), M(7,7);
    fullM.setZero();
    for( int i=0; i<m->nv; i++ ){
        for( int j=0; j<m->nv; j++ )
        if(j>=i)
            fullM(i,j) = denseM[i*m->nv + j];
    }
    M = fullM.block<7,7>(6,6);
    std::cout << "MuJoCo:\n" << M << "\n";

    // Calculate and print mass matrix differences
    std::cout << "Discrepancy:\n" << (M-data.M).cwiseAbs() << "\n";
    std::cout << "Max discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";

    // Shut down
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}