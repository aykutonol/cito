#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

// Perturbation flags
int pos_pert=1, vel_pert=0, tau_pert=0;

// Print flags
bool print_derivatives = false;

void copyData(const mjModel* m, const mjData* dmain, mjData* d)
{
    // copy state and control from dmain to thread-specific d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, m->nv);
    mju_copy(d->qacc, dmain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);
}

PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) getExternalForce(const mjModel* m, const mjData* d, const pinocchio::Model& model)
{
    PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext((size_t)model.njoints, pinocchio::Force::Zero());
    pinocchio::Force::Vector6 contact_force_ref;
    pinocchio::Model::JointIndex pin_jnt_id;
    Eigen::VectorXd vec_j2c_w(3), vec_j2c_b(3), fcon_w(3), fcon_b(3), hcon(6), hcon_neg(6), h_con_b_at_j(6);
    int geom_jntadr[2], geom_bodyid[2];
    for (int i=0; i<d->ncon; i++)
    {
        // Get contact force in contact frame
        mj_contactForce(m, d, i, hcon.data());
        hcon_neg = -hcon;
        // Get joint and body ids from the MuJoCo model
        geom_bodyid[0] = m->geom_bodyid[d->contact[i].geom1];
        geom_bodyid[1] = m->geom_bodyid[d->contact[i].geom2];
        geom_jntadr[0] = m->body_jntadr[geom_bodyid[0]];
        geom_jntadr[1] = m->body_jntadr[geom_bodyid[1]];
        // For each body involved in the contact, set the external force
        for(int j=0; j<2; j++)
        {
            if(geom_jntadr[j] != -1)
            {
                // Get the joint id from the Pinocchio model
                pin_jnt_id = model.getJointId(mj_id2name(m, mjOBJ_JOINT, geom_jntadr[j]));
                // Continue if the joint exists in the Pinocchio model
                if(pin_jnt_id < model.njoints)
                {
                    // Calculate vector from the joint anchor to the contact point
                    mju_sub3(vec_j2c_w.data(), d->contact[i].pos, d->xanchor+3*geom_jntadr[j]);
                    // Represent the contact wrench in the world frame
                    if(j==0)
                    {
                        mju_rotVecMatT(fcon_w.data(), hcon_neg.head(3).data(), d->contact[i].frame);
                    }
                    else
                    {
                        mju_rotVecMatT(fcon_w.data(), hcon.head(3).data(), d->contact[i].frame);
                    }
                    // Represent the contact wrench in the body frame at the contact point
                    mju_rotVecMatT(fcon_b.data(), fcon_w.head(3).data(), d->xmat+9*geom_bodyid[j]);
                    mju_rotVecMatT(vec_j2c_b.data(), vec_j2c_w.data(), d->xmat+9*geom_bodyid[j]);
                    h_con_b_at_j.head(3) = fcon_b;
                    h_con_b_at_j[3] = -vec_j2c_b[2]*fcon_b[1] + vec_j2c_b[1]*fcon_b[2];
                    h_con_b_at_j[4] =  vec_j2c_b[2]*fcon_b[0] - vec_j2c_b[0]*fcon_b[2];
                    h_con_b_at_j[5] = -vec_j2c_b[1]*fcon_b[0] + vec_j2c_b[0]*fcon_b[1];
                    // Set joint forces due to the contact
                    contact_force_ref = h_con_b_at_j;
                    fext[pin_jnt_id] += pinocchio::ForceRef<pinocchio::Force::Vector6>(contact_force_ref);
                }
            }
        }
    }
    return fext;
}

void printContactInfo(const mjModel* m, const mjData* d)
{
    Eigen::VectorXd hcon(6);
    std::cout << "\nContacts:\n";
    for(int i=0; i<d->ncon; i++)
    {
        mj_contactForce(m, d, i, hcon.data());
        int geom1_id=d->contact[i].geom1, geom2_id=d->contact[i].geom2;
        std::cout << "\tContact " << i << ", pos: " << d->contact[i].pos[0] << ", " <<
                        d->contact[i].pos[1] << ", " << d->contact[i].pos[2] << ", " <<
                        "dist: " << d->contact[i].dist << ", dim: " << d->contact[i].dim <<
                        "\n\t\tgeom1: " << mj_id2name(m, mjOBJ_GEOM, geom1_id) <<
                        ", geom2: " << mj_id2name(m, mjOBJ_GEOM, geom2_id) <<
                        ", hcon: " << hcon.head(3).transpose() << "\n";
        std::cout << "\t\tgeom1 body: " << mj_id2name(m, mjOBJ_BODY, m->geom_bodyid[geom1_id]) <<
                        ", jntnum: " << m->body_jntnum[m->geom_bodyid[geom1_id]] << 
                        ", jntadr: " << m->body_jntadr[m->geom_bodyid[geom1_id]] <<
                        ", dofnum: " << m->body_dofnum[m->geom_bodyid[geom1_id]] << 
                        ", dofadr: " << m->body_dofadr[m->geom_bodyid[geom1_id]] <<
                        "\n\t\tgeom2 body: " << mj_id2name(m, mjOBJ_BODY, m->geom_bodyid[geom2_id]) << 
                        ", jntnum: " << m->body_jntnum[m->geom_bodyid[geom2_id]] << 
                        ", jntadr: " << m->body_jntadr[m->geom_bodyid[geom2_id]] <<
                        ", dofnum: " << m->body_dofnum[m->geom_bodyid[geom2_id]] << 
                        ", dofadr: " << m->body_dofadr[m->geom_bodyid[geom2_id]] <<"\n";
    }
    if(d->ncon==0)
        std::cout << "\tNo contact\n";
}

int main(int argc, char const *argv[]) {
    // Get the perturbation flags, if set
    if(argc>1)
    {
        pos_pert = atoi(argv[1]);
        vel_pert = atoi(argv[2]);
        tau_pert = atoi(argv[3]);
    }
    // MuJoCo
    // Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Load model
    std::string mjModelPathStr = paths::workspaceDir + "/src/cito/model/sawyer_push.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    // Create data
    mjData* d = mj_makeData(m);
    
    // Set to a contact instant for Sawyer box pushing model
    // Position (from t=0.070 s)
    d->qpos[0] = 1.10264;
    d->qpos[1] = -0.000452574;
    d->qpos[2] = 0.0660566;
    d->qpos[3] = 0.999997;
    d->qpos[4] = -0.00127484;
    d->qpos[5] = -0.000904901;
    d->qpos[6] = 0.00206843;
    d->qpos[7] = -0.2432;
    d->qpos[8] = -0.5791;
    d->qpos[9] = -0.0017;
    d->qpos[10] = 1.6396;
    d->qpos[11] = -0.4296;
    d->qpos[12] = -0.4659;
    d->qpos[13] = 0.2362;
    // Velocity  (from t=0.070 s)
    d->qvel[0] = 0.528398;
    d->qvel[1] = -0.0905147;
    d->qvel[2] = 0.325857;
    d->qvel[3] = -0.509935;
    d->qvel[4] = -0.361961;
    d->qvel[5] = 0.827372;
    d->qvel[6] = 0.2546;
    d->qvel[7] = 0.9922;
    d->qvel[8] = 0.0203;
    d->qvel[9] = -4.4834;
    d->qvel[10] = -0.3477;
    d->qvel[11] = 7.0372;
    d->qvel[12] = 6.8713;
    // Control (from t=0.075 s)
    d->ctrl[0] = -2.04539;
    d->ctrl[1] = -26.0581;
    d->ctrl[2] = -4.33632;
    d->ctrl[3] = 4.06073;
    d->ctrl[4] = -1.64713;
    d->ctrl[5] = -2.80625;
    d->ctrl[6] = -0.100421;

    // Evaluate the forward dynamics, but do not integrate
    mj_forward(m, d);

    // Pinocchio
    // Create model & data
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);

    // Offsets between the models assuming the MuJoCo model may have more DOF
    const int ndof=model.nv, pos_off=m->nq-model.nq, vel_off=m->nv-model.nv;

    // Set the configuration identical to MuJoCo
    Eigen::VectorXd q(model.nq), v(model.nv), tau(model.nv), tau_w_contact(model.nv);
    mju_copy(q.data(), d->qpos+pos_off, model.nq);
    mju_copy(v.data(), d->qvel+vel_off, model.nv);
    mju_copy(tau.data(), d->ctrl, model.nv);
    // mju_copy(tau.data(), d->qfrc_applied+vel_off, model.nv);
    
    std::cout << "qacc after forward dynamics evaluation: ";
    mju_printMat(d->qacc, 1, m->nv);

    // Print configuration
    std::cout << "qpos:     ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel:     ";
    mju_printMat(d->qvel, 1, m->nv);
    std::cout << "ctrl:     ";
    mju_printMat(d->ctrl, 1, m->nu);
    std::cout << "qfrc:     ";
    mju_printMat(d->qfrc_applied, 1, m->nv);
    std::cout << "qext:     ";
    mju_printMat(d->qfrc_constraint, 1, m->nv);
    std::cout << "qfrc_act: ";
    mju_printMat(d->qfrc_actuator, 1, m->nv);

    // Print contact info
    printContactInfo(m, d);
    // Get contact forces from MuJoCo and set fext
    PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext;
    fext = getExternalForce(m, d, model);

    // Evaluate forward dynamics
    // MuJoCo
    mj_forward(m, d);
    // Pinocchio
    pinocchio::computeAllTerms(model, data, q, v);
    pinocchio::aba(model, data, q, v, tau, fext);

    // Get contact forces projected onto joint space
    Eigen::VectorXd qcon(ndof);
    mju_copy(qcon.data(), d->qfrc_constraint+vel_off, ndof);
    tau_w_contact = tau + qcon;
    
    // MuJoCo data
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
    M = fullM.block(vel_off,vel_off, ndof, ndof);
    // Bias term
    Eigen::VectorXd mj_bias(ndof);
    mju_copy(mj_bias.data(), d->qfrc_bias+vel_off, ndof);
    // Acceleration
    Eigen::VectorXd mj_qacc(ndof), mj_qacc_unc(ndof), pn_a(ndof), pn_a_unc(ndof);
    mju_copy(mj_qacc.data(), d->qacc+vel_off, ndof);
    mju_copy(mj_qacc_unc.data(), d->qacc_unc+vel_off, ndof);
    pn_a = data.ddq;

    // Compare mass matrices
    std::cout << "\nMass matrices:\n";
    // std::cout << "Pinocchio:\n" << data.M << "\n";
    // std::cout << "MuJoCo:\n" << M << "\n";
    // std::cout << "Discrepancy:\n" << (M-data.M).cwiseAbs() << "\n";
    std::cout << "Max discrepancy: " << ((M-data.M).cwiseAbs()).maxCoeff() << "\n";

    // Compare bias terms
    std::cout << "\nBias terms:\n";
    std::cout << "Pinocchio: " << data.nle.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_bias.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_bias-data.nle).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_bias-data.nle).cwiseAbs()).maxCoeff() << "\n";

    // Compare accelerations
    std::cout << "\nAccelerations:\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Unconstrained accelerations
    pinocchio::aba(model, data, q, v, tau);
    pn_a_unc = data.ddq;
    std::cout << "\nUnconstrained accelerations:\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc_unc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc_unc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc_unc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Derivative comparison
    // Create CITO class objects for MuJoCo calculations
    CitoParams cp(m);
    CitoNumDiff nd(m);
    // Initialize derivative matrices
    eigMd dqacc_dqpos, dqacc_dqvel, dqacc_dctrl,
          mj_da_dq, mj_da_dv, mj_da_du,
          pn_da_dq, pn_da_dv, pn_da_du,
          pn_da_dq_wo_fext, pn_da_dv_wo_fext, pn_da_du_wo_fext;
    dqacc_dqpos.resize(m->nv, m->nv);    dqacc_dqvel.resize(m->nv, m->nv);    dqacc_dctrl.resize(m->nv, m->nu); 
    mj_da_dq.resize(model.nv, model.nv); mj_da_dv.resize(model.nv, model.nv); mj_da_du.resize(model.nv, model.nv);
    pn_da_dq.resize(model.nv, model.nv); pn_da_dv.resize(model.nv, model.nv); pn_da_du.resize(model.nv, model.nv); 
    pn_da_dq_wo_fext.resize(model.nv, model.nv);
    pn_da_dv_wo_fext.resize(model.nv, model.nv);
    pn_da_du_wo_fext.resize(model.nv, model.nv); 
    
    // MuJoCo
    mjtNum* deriv = 0;
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);
    nd.worker(d, deriv);
    mju_copy(dqacc_dqpos.data(), deriv, m->nv*m->nv);
    mju_copy(dqacc_dqvel.data(), deriv+m->nv*m->nv, m->nv*m->nv);
    mju_copy(dqacc_dctrl.data(), deriv+2*m->nv*m->nv, m->nv*m->nu);
    mju_free(deriv);
    mj_da_dq = dqacc_dqpos.bottomRightCorner(model.nv, model.nv);
    mj_da_dv = dqacc_dqvel.bottomRightCorner(model.nv, model.nv);
    mj_da_du = dqacc_dctrl.bottomRows(model.nv);

    // Pinocchio w/ fext
    pinocchio::computeABADerivatives(model, data, q, v, tau, fext);
    pn_da_dq = data.ddq_dq;
    pn_da_dv = data.ddq_dv;
    pn_da_du = data.Minv;

    // Pinocchio w/o fext
    pinocchio::computeABADerivatives(model, data, q, v, tau_w_contact);
    pn_da_dq_wo_fext = data.ddq_dq;
    pn_da_dv_wo_fext = data.ddq_dv;
    pn_da_du_wo_fext = data.Minv;

    // Print derivative matrices
    if(print_derivatives)
        std::cout << "\nda_dq:\nMuJoCo:\n" << mj_da_dq << "\nPinocchio w/ fext:\n" << pn_da_dq <<
                    "\nPinocchio w/o fext:\n" << pn_da_dq_wo_fext <<
                    "\nda_dv:\nMuJoCo:\n" << mj_da_dv << "\nPinocchio w/ fext:\n" << pn_da_dv <<
                    "\nPinocchio w/o fext:\n" << pn_da_dv_wo_fext <<
                    "\nda_du:\nMuJoCo:\n" << mj_da_du << "\nPinocchio w/ fext:\n" << pn_da_du <<
                    "\nPinocchio w/o fext:\n" << pn_da_du_wo_fext << "\n";

    // Prediction accuracy comparison
    // Set random seed
    std::srand(std::time(0));
    // Generate random perturbation
    eigVd dq(ndof), dv(ndof), du(ndof);
    dq.setZero(); dv.setZero(); du.setZero();
    if(pos_pert)
        dq = Eigen::VectorXd::Random(ndof)*5e-2;
    if(vel_pert)
        dv = Eigen::VectorXd::Random(ndof)*5e-2;
    if(tau_pert)
        du = Eigen::VectorXd::Random(ndof)*5e-2;
    std::cout << "\nPerturbations:\n\tdq: " << dq.transpose() << "\n\tdv: " << 
                 dv.transpose() << "\n\tdu: " << du.transpose() << "\n";
    // Perturb states and controls
    q += dq;
    v += dv;
    tau += du;
    mjData* dPerturbed = mj_makeData(m);
    copyData(m, d, dPerturbed);
    mju_copy(dPerturbed->qpos+pos_off, q.data(), ndof);
    mju_copy(dPerturbed->qvel+vel_off, v.data(), ndof);
    mju_copy(dPerturbed->ctrl, tau.data(), ndof);

    // Calculate actual perturbed accelerations using MuJoCo
    Eigen::VectorXd mj_qacc_pert(ndof), mj_qacc_unc_pert(ndof), pn_a_pert(ndof), pn_a_unc_pert(ndof);
    mj_forward(m, dPerturbed);
    mju_copy(mj_qacc_pert.data(), dPerturbed->qacc+vel_off, ndof);
    mju_copy(mj_qacc_unc_pert.data(), dPerturbed->qacc_unc+vel_off, ndof);
    // Get constraint forces in joint space
    mju_copy(qcon.data(), dPerturbed->qfrc_constraint+vel_off, ndof);
    tau_w_contact = tau + qcon;

    // Print updated contact info
    printContactInfo(m, dPerturbed);
    // Update fext for Pinocchio computations
    fext = getExternalForce(m, dPerturbed, model);

    // Calculate actual perturbed accelerations using Pinocchio
    pinocchio::aba(model, data, q, v, tau, fext);
    pn_a_pert = data.ddq;
    pinocchio::aba(model, data, q, v, tau_w_contact);
    pn_a_unc_pert = data.ddq;

    // Calculate predicted accelerations
    Eigen::VectorXd mj_a_pred(ndof), pn_a_pred(ndof), pn_a_pred_wo_fext(ndof);
    mj_a_pred = mj_qacc + mj_da_dq*dq + mj_da_dv*dv + mj_da_du*du;
    pn_a_pred = mj_qacc + pn_da_dq*dq + pn_da_dv*dv + pn_da_du*du;
    pn_a_pred_wo_fext = mj_qacc + pn_da_dq_wo_fext*dq + pn_da_dv_wo_fext*dv + pn_da_du_wo_fext*du;

    // Print predicted accelerations
    std::cout << "\n\nActual accelerations after perturbation:" <<
                 "\nMuJoCo:             " << mj_qacc_pert.transpose() <<
                 "\nPinocchio:          " << pn_a_pert.transpose() <<
                 "\nPinocchio w/o fext: " << pn_a_unc_pert.transpose() <<
                 "\nPredicted accelerations:" <<
                 "\nMuJoCo:             " << mj_a_pred.transpose() <<
                 "\nPinocchio:          " << pn_a_pred.transpose() <<
                 "\nPinocchio w/o fext: " << pn_a_pred_wo_fext.transpose() <<
                 "\nDiscrepancies w.r.t. MuJoCo's perturbed acceleration:" <<
                 "\nMuJoCo:             " << ((mj_qacc_pert-mj_a_pred).cwiseAbs()).transpose() <<
                 "\nPinocchio:          " << ((mj_qacc_pert-pn_a_pred).cwiseAbs()).transpose() <<
                 "\nPinocchio w/o fext: " << ((mj_qacc_pert-pn_a_pred_wo_fext).cwiseAbs()).transpose() << "\n";

    // Shut down
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}