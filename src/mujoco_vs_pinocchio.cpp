#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

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
    // Set random velocity
    std::srand(std::time(0));
    Eigen::VectorXd vRand(m->nv);
    vRand.setZero();
    vRand.tail(7) = Eigen::VectorXd::Random(7);
    mju_copy(d->qvel, vRand.data(), m->nv);
    std::cout << "qpos: ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel: ";
    mju_printMat(d->qvel, 1, m->nv);
    /// Evaluate forward dynamics
    mj_forward(m, d);
    /// Get mass matrix
    mjtNum* denseM = mj_stackAlloc(d, m->nv*m->nv);
    mj_fullM(m, denseM, d->qM);
    Eigen::MatrixXd fullM(m->nv, m->nv), M(7, 7);
    fullM.setZero();
    for( int i=0; i<m->nv; i++ ){
        for( int j=0; j<m->nv; j++ )
        if(j>=i)
            fullM(i,j) = denseM[i*m->nv + j];
    }
    M = fullM.block<7,7>(6,6);
    /// Get bias term
    Eigen::VectorXd mj_bias(7);
    mju_copy(mj_bias.data(), d->qfrc_bias+6, 7);

    // Pinocchio
    /// Create model & data
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);
    /// Set a configuration identical to MuJoCo
    Eigen::VectorXd q(model.nq), v(model.nv);
    mju_copy(q.data(), d->qpos+7, model.nq);
    mju_copy(v.data(), d->qvel+6, model.nq);
    /// Evaluate forward dynamics
    pinocchio::computeAllTerms(model, data, q, v);

    // Print mass matrices
    std::cout << "Pinocchio:\n" << data.M << "\n";
    std::cout << "MuJoCo:\n" << M << "\n";
    std::cout << "Discrepancy:\n" << (M-data.M).cwiseAbs() << "\n";
    std::cout << "Max discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";

    // Compare bias terms
    std::cout << "\nBias terms:\nPinocchio:\t" << data.nle.transpose() << "\n";
    std::cout << "MuJoCo:\t" << mj_bias.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_bias-data.nle).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_bias-data.nle).cwiseAbs()).maxCoeff() << "\n";

    // Shut down
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}