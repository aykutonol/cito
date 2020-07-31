#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/model.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

// Perturbation flags
int pos_pert=1, vel_pert=0, tau_pert=0;

// Print flags
bool print_derivatives=true;

// Model flag
int include_object_in_pin=1;

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
        pos_pert = atoi(argv[1]);
    if(argc>2)
        vel_pert = atoi(argv[2]);
    if(argc>3)
        tau_pert = atoi(argv[3]);
    if(argc>4)
        include_object_in_pin = atoi(argv[4]);
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
    // d->qpos[2] = 0.064;      // for initial table contact
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
    // Free body translational velocities in global coordinates
    d->qvel[0] = 0.528398;
    d->qvel[1] = -0.0905147;
    d->qvel[2] = 0.325857;
    // Free body angular velocities in local coordinates
    d->qvel[3] = -0.509935;
    d->qvel[4] = -0.361961;
    d->qvel[5] = 0.827372;
    // Arm joint velocities
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

    // Evaluate the forward dynamics using MuJoCo
    mj_forward(m, d);

    // Pinocchio model & data
    pinocchio::Model model, model_obj, model_rbt;
    if(include_object_in_pin==0)
        pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model);
    else
    {
        pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/sawyer.urdf", model_rbt);
        pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/box.urdf", model_obj);
        // pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/box.urdf", pinocchio::JointModelFreeFlyer(), model_obj);
        model_obj.frames[1].name = "box_root_joint";
        pinocchio::appendModel(model_rbt, model_obj, 0, pinocchio::SE3::Identity(), model);
    }
    pinocchio::Data data(model);

    // Offsets between the models assuming the MuJoCo model may have initial unactuated DOF
    const int ndof=model.nv, pos_off=m->nq-model.nq, vel_off=m->nv-model.nv;

    // Assume a free joint is added if nq = nv+1 in Pinocchio
    bool has_free_jnt = (model.nq == model.nv + 1) ? true : false;
    std::cout << "\n\tINFO: Pinocchio model includes " << has_free_jnt << " free joint.\n\n";

    // Set the configuration identical to MuJoCo
    Eigen::VectorXd q(model.nq), v(model.nv), tau(model.nv), tau_w_contact(model.nv);
    mju_copy(q.data(), d->qpos+pos_off, model.nq);
    mju_copy(v.data(), d->qvel+vel_off, model.nv);
    mju_copy(tau.data(), d->qfrc_actuator+vel_off, model.nv);

    // Get object position and DOF indices from Pinocchio and convert MuJoCo free-body
    // states into Pinocchio's convention
    int obj_jnt_id=-1, obj_idx_v=-1;
    mjtNum q_g2l[4], q_l2g[4];
    if(has_free_jnt)
    {
        obj_jnt_id=model.getJointId("object");
        obj_idx_v=model.joints[obj_jnt_id].idx_v();
        for(int i=0; i<model.njoints; i++)
        {
            int idx_v=model.joints[i].idx_v(),
                idx_q=model.joints[i].idx_q();
            if(idx_v>=0)
            {
                // Check if the joint is free
                if(model.joints[i].nv()==6)
                {
                    // Change quaternion convention from wxyz to xyzw
                    for(int j=0; j<3; j++)
                    {
                        q[idx_q+j] = d->qpos[pos_off+idx_q+j];
                        q[idx_q+3+j] = d->qpos[pos_off+idx_q+4+j];
                    }
                    q[idx_q+6] = d->qpos[pos_off+idx_q+3];
                    // Project free-body linear velocity onto the local frame
                    Eigen::VectorXd vobj_g(3), vobj_l(3);
                    mju_copy3(vobj_g.data(), d->qvel+vel_off+idx_v);
                    mju_copy4(q_l2g, d->qpos+pos_off+idx_q+3);
                    mju_negQuat(q_g2l, q_l2g);
                    mju_rotVecQuat(vobj_l.data(), vobj_g.data(), q_g2l);
                    v.segment(idx_v, 3) = vobj_l;
                }
            }
        }
    }

    // Print MuJoCo data
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

    // Evaluate the forward dynamics using Pinocchio
    pinocchio::computeAllTerms(model, data, q, v);
    pinocchio::aba(model, data, q, v, tau, fext);

    // Print Pinocchio data
    std::cout << "Pinocchio:\nq: " << q.transpose() <<
                 "\nv: " << v.transpose() <<
                 "\na: " << data.ddq.transpose() <<
                 "\nu: " << tau.transpose() << "\n";

    // Get contact forces projected onto joint space
    Eigen::VectorXd qcon(ndof);
    mju_copy(qcon.data(), d->qfrc_constraint+vel_off, ndof);
    tau_w_contact = tau + qcon;

    // Project free-body bias/accelerations in Pinocchio onto the world frame
    Eigen::VectorXd pin_bias(model.nv), pn_a(model.nv), pn_a_unc(model.nv);
    pin_bias = data.nle;
    pn_a = data.ddq;
    // Replace linear spatial accelerations by classical
    pinocchio::Motion obj_acc;
    if(has_free_jnt)
    {
        pinocchio::forwardKinematics(model, data, q, v, data.ddq);
        obj_acc = pinocchio::getClassicalAcceleration(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        pn_a.segment(obj_idx_v, 3) = obj_acc.linear();
        // Compare spatial/local velocity and accelerations of the free body
        Eigen::VectorXd mj_obj_vel(6), mj_obj_acc(6), mj_obj_cvel(6), mj_obj_cacc(6);
        int mj_obj_bodyid = mj_name2id(m, mjOBJ_BODY, "object");
        mj_rnePostConstraint(m, d);
        mj_comVel(m, d);
        mj_objectVelocity(m, d, mjOBJ_BODY, mj_obj_bodyid, mj_obj_vel.data(), 0);
        mj_objectAcceleration(m, d, mjOBJ_BODY, mj_obj_bodyid, mj_obj_acc.data(), 0);
        mju_copy(mj_obj_cvel.data(), d->cvel+mj_obj_bodyid*6, 6);
        mju_copy(mj_obj_cacc.data(), d->cacc+mj_obj_bodyid*6, 6);
        pinocchio::Motion pn_obj_vel = pinocchio::getVelocity(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        pinocchio::Motion pn_obj_acc = pinocchio::getAcceleration(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        std::cout << "\nObject LWA velocity:\n\tPinocchio: w = " << pn_obj_vel.angular().transpose() << 
                    ", v =" << pn_obj_vel.linear().transpose() <<
                    "\n\tlocal:     w = " << mj_obj_vel.head(3).transpose() << 
                    ", v = " << mj_obj_vel.tail(3).transpose() <<
                    "\n\tcvel:     w = " << mj_obj_cvel.head(3).transpose() << 
                    ", v = " << mj_obj_cvel.tail(3).transpose() << "\n";
        std::cout << "\nObject LWA acceleration:\n\tPinocchio: w = " << pn_obj_acc.angular().transpose() << 
                    ", v =" << pn_obj_acc.linear().transpose() <<
                    "\n\tlocal:     w = " << mj_obj_acc.head(3).transpose() << 
                    ", v = " << mj_obj_acc.tail(3).transpose() <<
                    "\n\tcacc:      w = " << mj_obj_cacc.head(3).transpose() << 
                    ", v = " << mj_obj_cacc.tail(3).transpose() << "\n";
    }
    
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
    Eigen::VectorXd mj_qacc(ndof), mj_qacc_unc(ndof);
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
    std::cout << "Pinocchio: " << pin_bias.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_bias.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_bias-pin_bias).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_bias-pin_bias).cwiseAbs()).maxCoeff() << "\n";

    // Compare accelerations
    std::cout << "\nAccelerations:\nPinocchio:   " << pn_a.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc-pn_a).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc-pn_a).cwiseAbs()).maxCoeff() << "\n";

    // Unconstrained accelerations
    pinocchio::aba(model, data, q, v, tau);
    pn_a_unc = data.ddq;
    // Replace linear spatial accelerations by classical
    if(has_free_jnt)
    {
        pinocchio::forwardKinematics(model, data, q, v, data.ddq);
        obj_acc = pinocchio::getClassicalAcceleration(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        pn_a_unc.segment(obj_idx_v, 3) = obj_acc.linear();
    }
    std::cout << "\nUnconstrained accelerations:\nPinocchio: " << pn_a_unc.transpose() << "\n";
    std::cout << "MuJoCo:    " << mj_qacc_unc.transpose() << "\n";
    std::cout << "Discrepancy:\n" << (mj_qacc_unc-pn_a_unc).cwiseAbs().transpose() << "\n";
    std::cout << "Max discrepancy: " << ((mj_qacc_unc-pn_a_unc).cwiseAbs()).maxCoeff() << "\n";

    // Create CITO class objects for MuJoCo calculations
    CitoParams cp(m);
    CitoNumDiff nd(m);

    // Derivative comparison
    eigMd pn_da_dq, pn_da_dv, pn_da_du,
          pn_da_dq_wo_fext, pn_da_dv_wo_fext, pn_da_du_wo_fext;
    pn_da_dq.resize(model.nv, model.nv); pn_da_dv.resize(model.nv, model.nv); pn_da_du.resize(model.nv, model.nv); 
    pn_da_dq_wo_fext.resize(model.nv, model.nv);
    pn_da_dv_wo_fext.resize(model.nv, model.nv);
    pn_da_du_wo_fext.resize(model.nv, model.nv); 
    
    // MuJoCo
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dqacc_dqpos, dqacc_dqvel, dqacc_dctrl,
                                                                           mj_da_dq, mj_da_dv, mj_da_du;
    dqacc_dqpos.resize(m->nv, m->nv);    dqacc_dqvel.resize(m->nv, m->nv);    dqacc_dctrl.resize(m->nv, m->nv);
    mj_da_dq.resize(model.nv, model.nv); mj_da_dv.resize(model.nv, model.nv); mj_da_du.resize(model.nv, model.nv);
    mjtNum* deriv = 0;
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);
    nd.worker(d, deriv);
    mju_copy(dqacc_dqpos.data(), deriv, m->nv*m->nv);
    mju_copy(dqacc_dqvel.data(), deriv+m->nv*m->nv, m->nv*m->nv);
    mju_copy(dqacc_dctrl.data(), deriv+2*m->nv*m->nv, m->nv*m->nv);
    mju_free(deriv);
    mj_da_dq = dqacc_dqpos.bottomRightCorner(model.nv, model.nv);
    mj_da_dv = dqacc_dqvel.bottomRightCorner(model.nv, model.nv);
    mj_da_du = dqacc_dctrl.bottomRightCorner(model.nv, model.nv);

    // Pinocchio w/ fext
    pinocchio::computeABADerivatives(model, data, q, v, tau, fext);
    pn_da_dq = data.ddq_dq;
    pn_da_dv = data.ddq_dv;
    pn_da_du = data.Minv;
    if(has_free_jnt)
    {
        pn_da_dq.block(0, 0, 6, 6).setZero();
        pn_da_dv.block(0, 0, 6, 6).setZero();
        pn_da_dv.block(3, 3, 3, 3) = data.ddq_dv.block(obj_idx_v+3, obj_idx_v+3, 3, 3);
    }

    // Pinocchio w/o fext
    pinocchio::computeABADerivatives(model, data, q, v, tau);
    pn_da_dq_wo_fext = data.ddq_dq;
    pn_da_dv_wo_fext = data.ddq_dv;
    pn_da_du_wo_fext = data.Minv;
    if(has_free_jnt)
    {
        pn_da_dq_wo_fext.block(0, 0, 6, 6).setZero();
        pn_da_dv_wo_fext.block(0, 0, 6, 6).setZero();
        pn_da_dv_wo_fext.block(3, 3, 3, 3) = data.ddq_dv.block(obj_idx_v+3, obj_idx_v+3, 3, 3);
    }

    // Print derivative matrices
    std::cout << std::fixed;
    if(print_derivatives)
        // std::cout << "\nda_dq:\nMuJoCo:\n" << mj_da_dq.topLeftCorner(6,6) << "\nPinocchio w/ fext:\n" << pn_da_dq.topLeftCorner(6,6) <<
        //             "\nPinocchio w/o fext:\n" << pn_da_dq_wo_fext.topLeftCorner(6,6) <<
        //             "\nda_dv:\nMuJoCo:\n" << mj_da_dv.topLeftCorner(6,6) << "\nPinocchio w/ fext:\n" << pn_da_dv.topLeftCorner(6,6) << 
        //             "\nda_du:\nMuJoCo:\n" << mj_da_du.topLeftCorner(6,6) << "\nPinocchio w/ fext:\n" << pn_da_du.topLeftCorner(6,6) << "\n";
        std::cout << "\nda_dq:\nMuJoCo:\n" << mj_da_dq << "\nPinocchio w/ fext:\n" << pn_da_dq <<
                    "\nPinocchio w/o fext:\n" << pn_da_dq_wo_fext <<
                    "\nda_dv:\nMuJoCo:\n" << mj_da_dv << "\nPinocchio w/ fext:\n" << pn_da_dv << 
                    "\nda_du:\nMuJoCo:\n" << mj_da_du << "\nPinocchio w/ fext:\n" << pn_da_du << "\n";

    // Prediction accuracy comparison
    // Set random seed
    std::srand(std::time(0));
    auto dummy_rand = rand();
    // Generate random perturbation in global coordinates
    eigVd dq(ndof), dv(ndof), du(ndof);
    dq.setZero(); dv.setZero(); du.setZero();
    if(pos_pert)
        dq = Eigen::VectorXd::Random(ndof)*5e-2;
    if(vel_pert)
        dv = Eigen::VectorXd::Random(ndof)*5e-2;
    if(tau_pert)
        du.tail(m->nu) = Eigen::VectorXd::Random(m->nu)*5e-2;
    std::cout << "\nPerturbations:\n\tdq: " << dq.transpose() << "\n\tdv: " << 
                 dv.transpose() << "\n\tdu: " << du.transpose() << "\n";
    // Create perturbed data for MuJoCo
    mjData* dPerturbed = mj_makeData(m);
    copyData(m, d, dPerturbed);
    // Perturb positions
    if(has_free_jnt)
    {
        for(int i=vel_off; i<m->nv; i++)
        {
            // Get joint id for this dof
            int jid = m->dof_jntid[i];
            // Get quaternion address and dof position within quaternion (-1: not in quaternion)
            int quatadr = -1, dofpos = 0;
            if(m->jnt_type[jid]==mjJNT_FREE && i>=m->jnt_dofadr[jid]+3)
            {
                quatadr = m->jnt_qposadr[jid] + 3;
                dofpos = i - m->jnt_dofadr[jid] - 3;
            }
            // apply quaternion or simple perturbation
            if(quatadr>=0)
            {
                mjtNum angvel[3] = {0,0,0};
                angvel[dofpos] = dq[i-vel_off];
                mju_quatIntegrate(dPerturbed->qpos+pos_off+quatadr, angvel, 1);
            }
            else
                dPerturbed->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] += dq[i-vel_off];
        }
    }
    mju_copy(q.data(), dPerturbed->qpos+pos_off, model.nq);
    // Perturb velocities
    mju_copy(v.data(), d->qvel+vel_off, ndof);
    v += dv;
    mju_copy(dPerturbed->qvel+vel_off, v.data(), ndof);
    // Convert MuJoCo free-body states into Pinocchio's convention
    if(has_free_jnt)
    {
        for(int i=0; i<model.njoints; i++)
        {
            int idx_v=model.joints[i].idx_v(),
                idx_q=model.joints[i].idx_q();
            if(idx_v>=0)
            {
                // Check if the joint is free
                if(model.joints[i].nv()==6)
                {
                    // Change quaternion convention from wxyz to xyzw
                    for(int j=0; j<3; j++)
                    {
                        q[idx_q+j] = dPerturbed->qpos[pos_off+idx_q+j];
                        q[idx_q+3+j] = dPerturbed->qpos[pos_off+idx_q+4+j];
                    }
                    q[idx_q+6] = dPerturbed->qpos[pos_off+idx_q+3];
                    // Project free-body linear velocity onto the local frame
                    Eigen::VectorXd vobj_g(3), vobj_l(3);
                    mju_copy3(vobj_g.data(), dPerturbed->qvel+vel_off+idx_v);
                    mju_copy4(q_l2g, dPerturbed->qpos+pos_off+idx_q+3);
                    mju_negQuat(q_g2l, q_l2g);
                    mju_rotVecQuat(vobj_l.data(), vobj_g.data(), q_g2l);
                    v.segment(idx_v, 3) = vobj_l;
                }
            }
        }
    }
    // Perturb controls
    tau.tail(m->nu) += du.tail(m->nu);
    mju_copy(dPerturbed->ctrl, tau.tail(m->nu).data(), m->nu);
    // Evaluate forward dynamics
    for(int i=0; i<3; i++)
        mj_forward(m, dPerturbed);
    // Print changes
    Eigen::VectorXd q_n(m->nq), q_p(m->nq), q_d(m->nq);
    mju_copy(q_n.data(), d->qpos, m->nq);
    mju_copy(q_p.data(), dPerturbed->qpos, m->nq);
    q_d = q_p - q_n;
    Eigen::VectorXd v_n(m->nv), v_p(m->nv), v_d(m->nv);
    mju_copy(v_n.data(), d->qvel, m->nv);
    mju_copy(v_p.data(), dPerturbed->qvel, m->nv);
    v_d = v_p - v_n;
    Eigen::VectorXd u_n(m->nu), u_p(m->nu), u_d(m->nu);
    mju_copy(u_n.data(), d->ctrl, m->nu);
    mju_copy(u_p.data(), dPerturbed->ctrl, m->nu);
    u_d = u_p - u_n;
    std::cout << "\nqn:\n" << q_n.transpose() << 
                 "\ndq:\n" << dq.transpose() <<
                 "\nqp:\n" << q_p.transpose() <<
                 "\nqd:\n" << q_d.transpose() <<
                 "\nq: \n" << q.transpose() << "\n";
    std::cout << "\nvn:\n" << v_n.transpose() << 
                 "\ndv:\n" << dv.transpose() <<
                 "\nvp:\n" << v_p.transpose() <<
                 "\nvd:\n" << v_d.transpose() <<
                 "\nv :\n" << v.transpose() << "\n";
    std::cout << "\nun:\n" << u_n.transpose() << 
                 "\ndu:\n" << du.transpose() <<
                 "\nup:\n" << u_p.transpose() <<
                 "\nud:\n" << u_d.transpose() <<
                 "\nu :\n" << tau.transpose() << "\n";

    // Calculate actual perturbed accelerations using MuJoCo
    Eigen::VectorXd mj_qacc_pert(ndof), mj_qacc_unc_pert(ndof), pn_a_pert(ndof), pn_a_unc_pert(ndof);
    mju_copy(mj_qacc_pert.data(), dPerturbed->qacc+vel_off, ndof);
    mju_copy(mj_qacc_unc_pert.data(), dPerturbed->qacc_unc+vel_off, ndof);
    // Get constraint forces in joint space
    mju_copy(qcon.data(), dPerturbed->qfrc_constraint+vel_off, ndof);
    // Rotate the linear joint-space forces to local coordinates for the free joint
    if(has_free_jnt)
    {
        Eigen::Vector3d qcon_vobj;
        mju_rotVecQuat(qcon_vobj.data(), qcon.head(3).data(), q_g2l);
        qcon.head(3) = qcon_vobj;
    }
    tau_w_contact = tau + qcon;

    // Print updated contact info
    printContactInfo(m, dPerturbed);
    // Update fext for Pinocchio computations
    fext = getExternalForce(m, dPerturbed, model);

    // Calculate actual perturbed accelerations using Pinocchio
    pinocchio::aba(model, data, q, v, tau, fext);
    pn_a_pert = data.ddq;
    // Replace spatial accelerations with classical acceleration for the free joint
    if(has_free_jnt)
    {
        pinocchio::forwardKinematics(model, data, q, v, data.ddq);
        obj_acc = pinocchio::getClassicalAcceleration(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        pn_a_pert.segment(obj_idx_v, 3) = obj_acc.linear();
    }
    // Replace fext with external forces projected on the joint space
    pinocchio::aba(model, data, q, v, tau_w_contact);
    pn_a_unc_pert = data.ddq;
    // Replace spatial accelerations with classical acceleration for the free joint
    if(has_free_jnt)
    {
        pinocchio::forwardKinematics(model, data, q, v, data.ddq);
        obj_acc = pinocchio::getClassicalAcceleration(model, data, obj_jnt_id, pinocchio::LOCAL_WORLD_ALIGNED);
        pn_a_unc_pert.segment(obj_idx_v, 3) = obj_acc.linear();
    }

    // Calculate predicted accelerations
    Eigen::VectorXd mj_a_pred(ndof), pn_a_pred(ndof), pn_a_pred_wo_fext(ndof);
    mj_a_pred = mj_qacc + mj_da_dq*dq + mj_da_dv*dv + mj_da_du*du;
    pn_a_pred = mj_qacc + pn_da_dq*dq + pn_da_dv*dv + pn_da_du*du;
    pn_a_pred_wo_fext = mj_qacc + pn_da_dq_wo_fext*dq + pn_da_dv_wo_fext*dv + pn_da_du_wo_fext*du;

    // Print predicted accelerations
    std::cout << "\n\nActual accelerations after perturbation:" <<
                 "\nMuJoCo:             " << mj_qacc_pert.transpose() <<
                 "\nPinocchio:          " << pn_a_pert.transpose() <<
                 "\nPinocchio w/ qcon:  " << pn_a_unc_pert.transpose() <<
                 "\nPredicted accelerations:" <<
                 "\nMuJoCo:             " << mj_a_pred.transpose() <<
                 "\nPinocchio:          " << pn_a_pred.transpose() <<
                 "\nPinocchio w/o fext: " << pn_a_pred_wo_fext.transpose() <<
                 "\nDiscrepancies w.r.t. MuJoCo's perturbed acceleration:" <<
                 "\nMuJoCo:             " << ((mj_qacc_pert-mj_a_pred).cwiseAbs()).transpose() <<
                 "\nPinocchio:          " << ((mj_qacc_pert-pn_a_pred).cwiseAbs()).transpose() <<
                 "\nPinocchio w/o fext: " << ((mj_qacc_pert-pn_a_pred_wo_fext).cwiseAbs()).transpose() << "\n";

    // Delete perturbed data
    mj_deleteData(dPerturbed);

    // Shut down
    mj_deleteModel(m);
    mj_deleteData(d);
    mj_deactivate();
    return 0;
}