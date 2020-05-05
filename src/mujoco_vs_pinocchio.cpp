#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

const int pos_off=7, vel_off=6, ndof=7;

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
    // Set random position & velocity
    std::srand(std::time(0));
    std::cout << "Dummy random value: " << rand() << "\n";
    Eigen::VectorXd qRand(ndof), vRand(ndof), uRand(m->nv);
    qRand = Eigen::VectorXd::Random(ndof);
    vRand = Eigen::VectorXd::Random(ndof);
    uRand.setZero();
    uRand.tail(ndof) = Eigen::VectorXd::Random(ndof);
    mju_copy(d->qpos+pos_off, qRand.data(), ndof);
    mju_copy(d->qvel+vel_off, vRand.data(), ndof);
    // mju_copy(d->ctrl, uRand.data(), ndof);
    mju_add(d->qfrc_applied, d->qfrc_bias, uRand.data(), ndof);
    // Take a few steps
    for(int i=0; i<1; i++)
        mj_step(m, d);

    // Pinocchio
    /// Create model & data
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);
    /// Set the configuration identical to MuJoCo
    Eigen::VectorXd q(model.nq), v(model.nv), tau(model.nv);
    mju_copy(q.data(), d->qpos+pos_off, model.nq);
    mju_copy(v.data(), d->qvel+vel_off, model.nq);
    // mju_copy(tau.data(), d->ctrl, model.nq);
    mju_copy(tau.data(), d->qfrc_applied+vel_off, model.nq);

    // Print configuration
    std::cout << "qpos: ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel: ";
    mju_printMat(d->qvel, 1, m->nv);
    // std::cout << "ctrl: ";
    // mju_printMat(d->ctrl, 1, m->nu);
    std::cout << "qfrc: ";
    mju_printMat(d->qfrc_applied, 1, m->nv);
    std::cout << "qext: ";
    mju_printMat(d->qfrc_constraint, 1, m->nv);

    // Constraint info
    for(int i=0; i<d->nefc; i++) {
        std::cout << "Constraint " << i << ": " << d->efc_type[i] << "\n";
    }

    // Set joint forces due to contacts
    PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext((size_t)model.njoints, pinocchio::Force::Zero());
    pinocchio::Force::Vector3 tau_contact = pinocchio::Force::Vector3::Zero();
    for(int i=1; i<model.njoints; i++)
    {
        std::cout << "MuJoCo joint axis: " << m->jnt_axis[(i-1)*3] << m->jnt_axis[(i-1)*3+1] << m->jnt_axis[(i-1)*3+2] << "\n";
        for(int j=0; j<3; j++)
        {
            tau_contact[j] =  m->jnt_axis[(i-1)*3+j]*d->qfrc_constraint[vel_off+i-1];
        }
        fext[i].angular(tau_contact);
        std::cout << fext[i] << "\n";
    }

    // Evaluate forward dynamics
    // MuJoCo
    mj_forward(m, d);
    /// Pinocchio
    pinocchio::computeAllTerms(model, data, q, v);
    pinocchio::aba(model, data, q, v, tau, fext);

    // MuJoCo data
    /// Mass matrix
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
    /// Bias term
    Eigen::VectorXd mj_bias(ndof);
    mju_copy(mj_bias.data(), d->qfrc_bias+vel_off, ndof);
    /// Acceleration
    Eigen::VectorXd mj_qacc(ndof), mj_qacc_unc(ndof);;
    mju_copy(mj_qacc.data(), d->qacc+vel_off, ndof);
    mju_copy(mj_qacc_unc.data(), d->qacc_unc+vel_off, ndof);

    // Compare mass matrices
    std::cout << "\nMass matrices:\n";
    // std::cout << "Pinocchio:\n" << data.M << "\n";
    // std::cout << "MuJoCo:\n" << M << "\n";
    // std::cout << "Discrepancy:\n" << (M-data.M).cwiseAbs() << "\n";
    std::cout << "Max discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";

    // Compare bias terms
    std::cout << "\nBias terms:\n";
    // std::cout << "Pinocchio: " << data.nle.transpose() << "\n";
    // std::cout << "MuJoCo:    " << mj_bias.transpose() << "\n";
    // std::cout << "Discrepancy:\n" << (mj_bias-data.nle).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_bias-data.nle).cwiseAbs()).maxCoeff() << "\n";

    // Compare accelerations
    std::cout << "\nAccelerations:\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Unconstrained accelerations
    pinocchio::aba(model, data, q, v, tau);
    std::cout << "\nUnconstrained\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc_unc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc_unc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc_unc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Shut down
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}