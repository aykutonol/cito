#include "mujoco.h"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

const std::string workspaceDir = std::getenv("CITO_WS");
int nSample = 5;

int main(int argc, char const *argv[]) {
    // Parse command line arguments
    if(argc>1)
    {
        nSample = atoi(argv[1]);
    }

    // Set random seed
    std::srand(std::time(0));

    // Initialize MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Load model
    std::string mjModelPathStr = workspaceDir + "/src/cito/model/sawyer_contact.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");

    // Initialize Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel(workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);

    // Offsets between the models assuming the MuJoCo model may have more DOF
    const int ndof=model.nv, pos_off=m->nq-model.nq, vel_off=m->nv-model.nv;

    // Initialize state and control vectors
    Eigen::VectorXd q(model.nq), v(model.nv), tau(model.nv), qcon(model.nv);
    q.setZero(); v.setZero(); tau.setZero(); qcon.setZero();
    // Random configurations
    Eigen::VectorXd qRand(m->nq), vRand(m->nv), uRand(m->nu);
    qRand.setZero(); vRand.setZero(); uRand.setZero();

    // Test loop
    for( int k=0; k<nSample; k++ )
    {
        std::cout << "\n====================== Sample no: " << k+1 << " ======================\n";
        // Generate random configuration
        if( k>0 )
        {
            qRand = Eigen::VectorXd::Random(m->nq)*2;
            vRand = Eigen::VectorXd::Random(m->nv)*1;
            uRand = Eigen::VectorXd::Random(m->nu)*1;
        }

        // Create MuJoCo data
        mjData* d = mj_makeData(m);
        mju_copy(d->qpos, qRand.data(), m->nq);
        mju_copy(d->qvel, vRand.data(), m->nv);
        mju_copy(d->ctrl, uRand.data(), m->nu);
        // Create Pinocchio data
        mju_copy(q.data(),   d->qpos+pos_off, model.nq);
        mju_copy(v.data(),   d->qvel+vel_off, model.nv);
        mju_copy(tau.data(), d->ctrl, model.nv);
        
        // Run MuJoCo forward dynamics
        mj_forward(m, d);
        // Run Pinocchio forward dynamics
        pinocchio::computeAllTerms(model, data, q, v);
        pinocchio::aba(model, data, q, v, tau);

        // Contact force in joint space
        mju_copy(qcon.data(), d->qfrc_constraint+vel_off, ndof);
        std::cout << "qcon: " << qcon.transpose() << "\n";

        // Mass matrix
        mjtNum* denseM = mj_stackAlloc(d, m->nv*m->nv);
        mj_fullM(m, denseM, d->qM);
        Eigen::MatrixXd fullM(m->nv, m->nv), M(ndof, ndof);
        fullM.setZero();
        for( int i=0; i<m->nv; i++ ){
            for( int j=0; j<m->nv; j++ )
            if(j>=i)
                fullM(i,j) = denseM[i*m->nv + j];
        }
        M = fullM.block(vel_off, vel_off, ndof, ndof);

        // Bias term
        Eigen::VectorXd mj_bias(ndof);
        mju_copy(mj_bias.data(), d->qfrc_bias+vel_off, ndof);

        // Acceleration
        Eigen::VectorXd mj_qacc(ndof), mj_qacc_unc(ndof);
        mju_copy(mj_qacc.data(), d->qacc+vel_off, ndof);
        mju_copy(mj_qacc_unc.data(), d->qacc_unc+vel_off, ndof);
        
        // Compare mass matrices
        std::cout << "\nMass matrices:\nPinocchio:\n" << data.M <<
                     "\nMuJoCo:\n" << M <<
                     "\nDiscrepancy:\n" << (M-data.M).cwiseAbs() <<
                     "\nMax discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";
        // Compare bias terms
        std::cout << "\nBias terms:  " <<
                     "\nPinocchio:   " << data.nle.transpose() <<
                     "\nMuJoCo:      " << mj_bias.transpose() <<
                     "\nDiscrepancy: " << (mj_bias-data.nle).cwiseAbs().transpose() <<
                     "\nMax discrepancy: " << ((mj_bias-data.nle).cwiseAbs()).maxCoeff() << "\n";
        // Compare unconstrained accelerations
        std::cout << "\nUnconstrained accelerations:" <<
                     "\nPinocchio:   " << data.ddq.transpose() <<
                     "\nMuJoCo:      " << mj_qacc_unc.transpose() <<
                     "\nDiscrepancy: " << (mj_qacc_unc-data.ddq).cwiseAbs().transpose() <<
                     "\nMax discrepancy: " << ((mj_qacc_unc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

        // Evaluate and compare constrained accelerations
        pinocchio::aba(model, data, q, v, tau+qcon);
        std::cout << "\n Constrained accelerations:" <<
                     "\nPinocchio:   " << data.ddq.transpose() <<
                     "\nMuJoCo:      " << mj_qacc.transpose() <<
                     "\nDiscrepancy: " << (mj_qacc-data.ddq).cwiseAbs().transpose() <<
                     "\nMax discrepancy: " << ((mj_qacc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

        // Delete data
        mj_deleteData(d);
    }

    // Shut down MuJoCo
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}