#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

const int pos_off=0, vel_off=0, ndof=7;

int main(int argc, char const *argv[]) {
    // MuJoCo
    /// Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    /// Load model
    std::string mjModelPathStr = paths::workspaceDir + "/src/cito/model/sawyer.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    /// Create data
    mjData* d = mj_makeData(m);
    /// Set the joint positions to the preset configuration in the model
    mju_copy(d->qpos, m->key_qpos, m->nq);
    // Set random position & velocity
    std::srand(std::time(0));
    std::cout << "Dummy random value: " << rand() << "\n";
    Eigen::VectorXd qRand(ndof), vRand(ndof), uRand(ndof);
    qRand = Eigen::VectorXd::Random(ndof);
    vRand = Eigen::VectorXd::Random(ndof)*10;
    uRand = Eigen::VectorXd::Random(ndof)*10;
    mju_copy(d->qpos+pos_off, qRand.data(), ndof);
    mju_copy(d->qvel+vel_off, vRand.data(), ndof);
    mju_copy(d->qfrc_applied+vel_off, uRand.data(), ndof);
    std::cout << "qpos: ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel: ";
    mju_printMat(d->qvel, 1, m->nv);
    std::cout << "qfrc: ";
    mju_printMat(d->qfrc_applied, 1, m->nv);
    /// Evaluate forward dynamics
    mj_forward(m, d);
    /// Get mass matrix
    mjtNum* denseM = mj_stackAlloc(d, m->nv*m->nv);
    mj_fullM(m, denseM, d->qM);
    Eigen::MatrixXd fullM(m->nv, m->nv), M(ndof, ndof);
    fullM.setZero();
    for( int i=0; i<m->nv; i++ ){
        for( int j=0; j<m->nv; j++ )
        if(j>=i)
            fullM(i,j) = denseM[i*m->nv + j];
    }
    M = fullM.block<ndof,ndof>(vel_off,vel_off);
    /// Get bias term
    Eigen::VectorXd mj_bias(ndof);
    mju_copy(mj_bias.data(), d->qfrc_bias+vel_off, ndof);
    /// Get acceleration
    Eigen::VectorXd mj_qacc(ndof);
    mju_copy(mj_qacc.data(), d->qacc+vel_off, ndof);

    // Pinocchio
    /// Create model & data
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);
    /// Set the configuration identical to MuJoCo
    Eigen::VectorXd q(model.nq), v(model.nv), tau(model.nv);
    mju_copy(q.data(), d->qpos+pos_off, model.nq);
    mju_copy(v.data(), d->qvel+vel_off, model.nq);
    mju_copy(tau.data(), d->qfrc_applied+vel_off, model.nq);
    // tau.setZero();
    /// Evaluate forward dynamics
    pinocchio::computeAllTerms(model, data, q, v);
    pinocchio::aba(model, data, q, v, tau);

    // Compare mass matrices
    std::cout << "Pinocchio:\n" << data.M << "\n";
    std::cout << "MuJoCo:\n" << M << "\n";
    std::cout << "Discrepancy:\n" << (M-data.M).cwiseAbs() << "\n";
    std::cout << "Max discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";

    // Compare bias terms
    std::cout << "\nBias terms:\nPinocchio: " << data.nle.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_bias.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_bias-data.nle).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_bias-data.nle).cwiseAbs()).maxCoeff() << "\n";

    // Compare accelerations
    std::cout << "\nAccelerations:\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Shut down
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}