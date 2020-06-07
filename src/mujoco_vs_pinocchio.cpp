#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

const int pos_off=7, vel_off=6, ndof=7;
int pos_pert=1, vel_pert=1, acc_pert=1;

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

int main(int argc, char const *argv[]) {
    // Get the perturbation flags, if set
    if(argc>1)
    {
        pos_pert = atoi(argv[1]);
        vel_pert = atoi(argv[2]);
        acc_pert = atoi(argv[3]);
    }
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
    // /// Set the joint positions to the preset configuration in the model
    // mju_copy(d->qpos, m->key_qpos, m->nq);
    // // Set random position & velocity
    // std::srand(std::time(0));
    // std::cout << "Dummy random value: " << rand() << "\n";
    // Eigen::VectorXd qRand(ndof), vRand(ndof), uRand(m->nv);
    // qRand = Eigen::VectorXd::Random(ndof);
    // vRand = Eigen::VectorXd::Random(ndof);
    // uRand.setZero();
    // uRand.tail(ndof) = Eigen::VectorXd::Random(ndof);
    // mju_copy(d->qpos+pos_off, qRand.data(), ndof);
    // mju_copy(d->qvel+vel_off, vRand.data(), ndof);
    // // mju_copy(d->ctrl, uRand.data(), ndof);
    // mju_add(d->qfrc_applied, d->qfrc_bias, uRand.data(), ndof);
    // Set to a contact instant for Sawyer box pushing model
    /// Position (from t=0.070 s)
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
    /// Velocity  (from t=0.070 s)
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
    /// Control (from t=0.075 s)
    d->ctrl[0] = -2.04539;
    d->ctrl[1] = -26.0581;
    d->ctrl[2] = -4.33632;
    d->ctrl[3] = 4.06073;
    d->ctrl[4] = -1.64713;
    d->ctrl[5] = -2.80625;
    d->ctrl[6] = -0.100421;
    // Control ( (from t=0.070 s))
    // d->ctrl[0] = -0.0177;
    // d->ctrl[1] = -29.0172;
    // d->ctrl[2] = -4.4444;
    // d->ctrl[3] = 5.8086;
    // d->ctrl[4] = -1.9072;
    // d->ctrl[5] = -2.911;
    // d->ctrl[6] = -0.0996;
    /// Evaluate the dynamics, but do not integrate
    mj_forward(m, d);
    // Take a few steps
    // for(int i=0; i<3; i++)
        // mj_step(m, d);

    // Pinocchio
    /// Create model & data
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    pinocchio::Data data(model);
    /// Set the configuration identical to MuJoCo
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

    // Create position and force vectors
    Eigen::VectorXd pos_con(3), dir_con(3), pos_jnt(3), pos_bdy(3),
                    vec_j2c(3), vec_b2j(3), vec_b2c(3),
                    f_efc(3), f_con_w(3), 
                    h_con(6), h_con_w_at_j(6), h_con_w_at_b(6);
    // Constraint info
    std::cout << "Constraints:\n";
    for(int i=0; i<1; i++) {
        std::cout << "\tConstraint " << i << ", type: " << d->efc_type[i] <<
                     ", dist: " << d->contact[i].dist <<
                     ", geom1: " << mj_id2name(m, mjOBJ_GEOM, d->contact[i].geom1) <<
                     ", geom2: " << mj_id2name(m, mjOBJ_GEOM, d->contact[i].geom2) <<
                     "\n\tefc_force: " << d->efc_force[i] <<
                     ", efc_state: " << d->efc_state[i] <<
                     ", efc_id: " << d->efc_id[i] << "\n";
        // Get positions of contact, joint, and body frames
        for(int j=0; j<3; j++)
        {
            pos_con(j) = d->contact[i].pos[j];
            dir_con(j) = d->contact[i].frame[j];
            pos_jnt(j) = d->xanchor[mj_name2id(m, mjOBJ_JOINT, "right_j6")*3+j];
            pos_bdy(j) = d->xipos[mj_name2id(m, mjOBJ_BODY, "right_l6")*3+j];
        }
        f_efc = d->efc_force[i]*dir_con;
        // Calculate vectors between frames
        vec_j2c = pos_con - pos_jnt;
        vec_b2j = pos_jnt - pos_bdy;
        vec_b2c = pos_con - pos_bdy;
        // Get the contact wrench in the contact frame
        mj_contactForce(m, d, i, h_con.data());
        // Represent the contact wrench in the world frame
        mju_rotVecMatT(f_con_w.data(), h_con.head(3).data(), d->contact[i].frame);
        h_con_w_at_j.head(3) = f_con_w;
        h_con_w_at_j(3) = -vec_j2c[2]*f_con_w[1] + vec_j2c[1]*f_con_w[2];
        h_con_w_at_j(4) =  vec_j2c[2]*f_con_w[0] - vec_j2c[0]*f_con_w[2];
        h_con_w_at_j(5) = -vec_j2c[1]*f_con_w[0] + vec_j2c[0]*f_con_w[1];
        h_con_w_at_b.head(3) = f_con_w;
        h_con_w_at_b(3) = -vec_b2c[2]*f_con_w[1] + vec_b2c[1]*f_con_w[2];
        h_con_w_at_b(4) =  vec_b2c[2]*f_con_w[0] - vec_b2c[0]*f_con_w[2];
        h_con_w_at_b(5) = -vec_b2c[1]*f_con_w[0] + vec_b2c[0]*f_con_w[1];
        // Print positions of and vectors between frames and forces
        std::cout << "\tpos_con: " << pos_con.transpose() << 
                     ", pos_jnt: " << pos_jnt.transpose() <<
                     ", pos_bdy: " << pos_bdy.transpose() <<
                     "\n\tdir_con: " << dir_con.transpose() <<
                     ", vec_j2c: " << vec_j2c.transpose() <<
                     "\n\tfefc: " << f_efc.transpose() <<
                     ", fcon_w: " << f_con_w.transpose() <<
                     "\n\thcon: " << h_con.transpose() <<
                     "\n\th_con_w_at_j: " << h_con_w_at_j.transpose() <<
                     "\n\th_con_w_at_j: " << h_con_w_at_b.transpose() << "\n\n";
        // Force applied at joint frame
        Eigen::Matrix<double, 3, vel_off+ndof, Eigen::RowMajor> Jct, Jcr;
        mj_jac(m, d, Jct.data(), Jcr.data(), pos_jnt.data(), mj_name2id(m, mjOBJ_BODY, "right_l6"));
        std::cout << "Joint Jacobian\n" <<
                     "Jct:\n" << Jct.block<3,ndof>(0,vel_off) << "\nJcr:\n" << Jcr.block<3,ndof>(0,vel_off) <<
                     "\nqcon: " << 
                     (Jct.block<3,ndof>(0,vel_off).transpose()*h_con_w_at_j.head(3)+Jcr.block<3,ndof>(0,vel_off).transpose()*h_con_w_at_j.tail(3)).transpose() 
                     << "\n\n";
        // Force applied at body frame
        mj_jac(m, d, Jct.data(), Jcr.data(), pos_bdy.data(), mj_name2id(m, mjOBJ_BODY, "right_l6"));
        std::cout << "Body Jacobian\n" <<
                     "Jct:\n" << Jct.block<3,ndof>(0,vel_off) << "\nJcr:\n" << Jcr.block<3,ndof>(0,vel_off) <<
                     "\nqcon: " << 
                     (Jct.block<3,ndof>(0,vel_off).transpose()*h_con_w_at_b.head(3)+Jcr.block<3,ndof>(0,vel_off).transpose()*h_con_w_at_b.tail(3)).transpose() 
                     << "\n\n";
        // Force applied at contact frame
        mj_jac(m, d, Jct.data(), Jcr.data(), pos_con.data(), mj_name2id(m, mjOBJ_BODY, "right_l6"));
        std::cout << "Constraint Jacobian\n" <<
                     "Jct:\n" << Jct.block<3,ndof>(0,vel_off) << "\nJcr:\n" << Jcr.block<3,ndof>(0,vel_off) <<
                     "\nqcon: " << 
                     (Jct.block<3,ndof>(0,vel_off).transpose()*f_con_w).transpose() << "\n\n";
    }
    mju_printMat(d->efc_J, d->nefc, m->nv);
    Eigen::VectorXd qefc(m->nv);
    mj_mulJacTVec(m, d, qefc.data(), d->efc_force);
    std::cout << "qefc: " << qefc.transpose() << "\n\n";

    // Represent the contact wrench in the body frame at the contact point
    Eigen::VectorXd f_con_b(3), h_con_b_at_j(6), vec_b_j2c(3);
    mju_rotVecMatT(f_con_b.data(), f_con_w.head(3).data(), d->xmat+9*mj_name2id(m, mjOBJ_BODY, "right_l6"));
    mju_rotVecMatT(vec_b_j2c.data(), vec_j2c.data(), d->xmat+9*mj_name2id(m, mjOBJ_BODY, "right_l6"));
    h_con_b_at_j.head(3) = f_con_b;
    h_con_b_at_j(3) = -vec_b_j2c[2]*f_con_b[1] + vec_b_j2c[1]*f_con_b[2];
    h_con_b_at_j(4) =  vec_b_j2c[2]*f_con_b[0] - vec_b_j2c[0]*f_con_b[2];
    h_con_b_at_j(5) = -vec_b_j2c[1]*f_con_b[0] + vec_b_j2c[0]*f_con_b[1];

    // Compare the actual joint axis with the initial axis rotated by the orientation of the body frame
    Eigen::VectorXd jnt_axis(3), jnt_axis_der(3);
    mju_copy3(jnt_axis.data(), d->xaxis+3*mj_name2id(m, mjOBJ_JOINT, "right_j6"));
    mju_printMat(m->jnt_axis+3*mj_name2id(m, mjOBJ_JOINT, "right_j6"), 1, 3);
    mju_rotVecMat(jnt_axis_der.data(), m->jnt_axis+3*mj_name2id(m, mjOBJ_JOINT, "right_j6"),
                                       d->xmat+9*mj_name2id(m, mjOBJ_BODY, "right_l6"));
    std::cout << "\nJoint axis: " << jnt_axis.transpose() << "\nDerived: " << jnt_axis_der.transpose() << "\n\n";

    // Set joint forces due to contacts
    PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext((size_t)model.njoints, pinocchio::Force::Zero());
    pinocchio::Force::Vector6 contact_force_ref = h_con_b_at_j;
    const pinocchio::Model::JointIndex con_jnt_id = model.getJointId("right_j6");
    fext[con_jnt_id] = pinocchio::ForceRef<pinocchio::Force::Vector6>(contact_force_ref);
    std::cout << "fext:\n";
    for(int i=0; i<model.njoints; i++)
        std::cout << fext[i] << "\n";

    // Evaluate forward dynamics
    // MuJoCo
    mj_forward(m, d);
    /// Pinocchio
    pinocchio::computeAllTerms(model, data, q, v);
    pinocchio::aba(model, data, q, v, tau, fext);
    // tau_w_contact = tau + qefc.tail(model.nv);
    // pinocchio::aba(model, data, q, v, tau_w_contact);
    // std::cout << "tau w/ contact forces: " << tau_w_contact.transpose() << "\n";
    
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
    std::cout << "\nUnconstrained\nPinocchio: " << data.ddq.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc_unc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc_unc-data.ddq).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc_unc-data.ddq).cwiseAbs()).maxCoeff() << "\n";

    // Derivative comparison
    /// Create CITO class objects for MuJoCo calculations
    CitoParams cp(m);
    CitoNumDiff nd(m);
    /// Initialize derivative matrices
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
    pinocchio::computeABADerivatives(model, data, q, v, tau);
    pn_da_dq_wo_fext = data.ddq_dq;
    pn_da_dv_wo_fext = data.ddq_dv;
    pn_da_du_wo_fext = data.Minv;

    // Print derivative matrices
    std::cout << "\nda_dq:\nMuJoCo:\n" << mj_da_dq << "\nPinocchio w/ fext:\n" << pn_da_dq <<
                 "\nPinocchio w/o fext:\n" << pn_da_dq_wo_fext <<
                 "\nda_dv:\nMuJoCo:\n" << mj_da_dv << "\nPinocchio w/ fext:\n" << pn_da_dv <<
                 "\nPinocchio w/o fext:\n" << pn_da_dv_wo_fext <<
                 "\nda_du:\nMuJoCo:\n" << mj_da_du << "\nPinocchio w/ fext:\n" << pn_da_du <<
                 "\nPinocchio w/o fext:\n" << pn_da_du_wo_fext << "\n";

    
    // Prediction accuracy comparison
    /// Set random seed
    std::srand(std::time(0));
    /// Generate random perturbation
    eigVd dq(ndof), dv(ndof), du(ndof);
    dq.setZero(); dv.setZero(); du.setZero();
    if(pos_pert)
        dq = Eigen::VectorXd::Random(ndof)*5e-2;
    if(vel_pert)
        dv = Eigen::VectorXd::Random(ndof)*5e-2;
    if(acc_pert)
        du = Eigen::VectorXd::Random(ndof)*5e-2;
    std::cout << "Perturbations:\n\tdq: " << dq.transpose() << "\n\tdv: " << 
                 dv.transpose() << "\n\tdu: " << du.transpose() << "\n";
    /// Perturb states and controls
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

    // Update fext for Pinocchio forward kinematics calculation
    std::cout << "\nConstraints after perturbation:\n";
    for(int i=0; i<dPerturbed->nefc; i++)
    {
        std::cout << "\tConstraint " << i << ", type: " << dPerturbed->efc_type[i] <<
                    ", dist: " << dPerturbed->contact[i].dist <<
                    ", geom1: " << mj_id2name(m, mjOBJ_GEOM, dPerturbed->contact[i].geom1) <<
                    ", geom2: " << mj_id2name(m, mjOBJ_GEOM, dPerturbed->contact[i].geom2) <<
                    "\n\tefc_force: " << dPerturbed->efc_force[i] <<
                    ", efc_state: " << dPerturbed->efc_state[i] <<
                    ", efc_id: " << dPerturbed->efc_id[i] << "\n";
    }
    /// Get positions of contact, joint, and body frames
    for(int j=0; j<3; j++)
    {
        pos_con(j) = dPerturbed->contact[0].pos[j];
        pos_jnt(j) = dPerturbed->xanchor[mj_name2id(m, mjOBJ_JOINT, "right_j6")*3+j];
    }
    /// Calculate vector from the joint anchor to the contact point
    vec_j2c = pos_con - pos_jnt;
    /// Get the contact wrench in the contact frame
    mj_contactForce(m, dPerturbed, 0, h_con.data());
    std::cout << "\nhcon: " << h_con.transpose() << "\n";
    /// Represent the contact wrench in the world frame
    mju_rotVecMatT(f_con_w.data(), h_con.head(3).data(), dPerturbed->contact[0].frame);
    /// Represent the contact wrench in the body frame at the contact point
    mju_rotVecMatT(f_con_b.data(), f_con_w.head(3).data(), dPerturbed->xmat+9*mj_name2id(m, mjOBJ_BODY, "right_l6"));
    mju_rotVecMatT(vec_b_j2c.data(), vec_j2c.data(), dPerturbed->xmat+9*mj_name2id(m, mjOBJ_BODY, "right_l6"));
    h_con_b_at_j.head(3) = f_con_b;
    h_con_b_at_j(3) = -vec_b_j2c[2]*f_con_b[1] + vec_b_j2c[1]*f_con_b[2];
    h_con_b_at_j(4) =  vec_b_j2c[2]*f_con_b[0] - vec_b_j2c[0]*f_con_b[2];
    h_con_b_at_j(5) = -vec_b_j2c[1]*f_con_b[0] + vec_b_j2c[0]*f_con_b[1];
    /// Set joint forces due to contacts
    contact_force_ref = h_con_b_at_j;
    fext[con_jnt_id] = pinocchio::ForceRef<pinocchio::Force::Vector6>(contact_force_ref);

    // Calculate actual perturbed accelerations using Pinocchio
    pinocchio::aba(model, data, q, v, tau, fext);
    pn_a_pert = data.ddq;
    pinocchio::aba(model, data, q, v, tau);
    pn_a_unc_pert = data.ddq;

    // Calculate predicted accelerations
    Eigen::VectorXd mj_a_pred(ndof), pn_a_pred(ndof), pn_a_pred_wo_fext(ndof);
    mj_a_pred = mj_qacc + mj_da_dq*dq + mj_da_dv*dv + mj_da_du*du;
    pn_a_pred = mj_qacc + pn_da_dq*dq + pn_da_dv*dv + pn_da_du*du;
    pn_a_pred_wo_fext = mj_qacc + pn_da_dq_wo_fext*dq + pn_da_dv_wo_fext*dv + pn_da_du_wo_fext*du;

    /// Print predicted accelerations
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