#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/compute-all-terms.hpp"

// Perturbation flags
int pos_pert=1, vel_pert=0, tau_pert=0;

// Compensate bias term
double compensateBias = 1.0;

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

void printConfig(mjModel* m, mjData* d)
{
    std::cout << "\ntime: " << d->time << "\n";
    std::cout << "qpos: ";
    mju_printMat(d->qpos, 1, m->nq);
    std::cout << "qvel: ";
    mju_printMat(d->qvel, 1, m->nv);
    std::cout << "ctrl: ";
    mju_printMat(d->ctrl, 1, m->nu);
    std::cout << "qfrc: ";
    mju_printMat(d->qfrc_applied, 1, m->nv);
    std::cout << "qacc: ";
    mju_printMat(d->qacc, 1, m->nv);
    std::cout << "qcst: ";
    mju_printMat(d->qfrc_constraint, 1, m->nv);
}

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

// Flags
bool printDeriv = false;
bool printMInit = false;
bool printPInit = false;
bool printConst = false;
bool printPred  = true;
bool printPert  = false;
bool printTraj  = false;
bool printTime  = false;

int main(int argc, char const *argv[]) {
    // Parse command line arguments
    int nSample = 50, fext_flag=1;
    if(argc>1)
    {
        nSample = atoi(argv[1]);
        if(argc>2)
        {
            fext_flag = atoi(argv[2]);
            // Get the perturbation flags, if set
            if(argc>3)
            {
                pos_pert = atoi(argv[3]);
                vel_pert = atoi(argv[4]);
                tau_pert = atoi(argv[5]);
            }
        }
    }

    // Initialize MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Load model
    std::string mjModelPathStr = paths::workspaceDir + "/src/cito/model/ur3e_contact.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    // Allocate derivatives
    mjtNum* deriv = 0;
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);

    // Assert all DOF are actuated
    assert(m->nv==m->nu);

    // Initialize Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/ur3e.urdf", model);
    pinocchio::Data data(model);

    // Assert the models have identical DOF
    assert(m->nq==model.nq);
    assert(m->nv==model.nv);
    assert(m->nu==model.nv);

    // Initialize CITO objects
    CitoParams  cp(m);
    CitoNumDiff nd(m);
    CitoControl cc(m);

    // MuJoCo state and control vectors
    eigVm x, u;
    x.resize(cp.n); u.resize(cp.m);
    x.setZero();    u.setZero();
    // Pinocchio state and control vectors
    eigVd q, v, tau, qcon;
    q.resize(model.nq); v.resize(model.nv); tau.resize(model.nv); qcon.resize(model.nv);
    q.setZero();        v.setZero();        tau.setZero();        qcon.setZero();

    // Derivative matrices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> da_dq, da_dv, da_df;
    da_dq.resize(m->nv, m->nv); da_dv.resize(m->nv, m->nv); da_df.resize(m->nv, m->nv);
    eigMd FxHW, FuHW, FxW, FuW, FxP, FuP;
    FxHW.resize(cp.n, cp.n);    FxW.resize(cp.n, cp.n);     FxP.resize(cp.n, cp.n);
    FuHW.resize(cp.n, cp.m);    FuW.resize(cp.n, cp.m);     FuP.resize(cp.n, cp.m);
    // Random configurations
    eigVd qRand(model.nq), vRand(model.nv), uRand(model.nv);
    qRand.setZero(); vRand.setZero(); uRand.setZero();
    // Perturbations
    eigVm dx(cp.n), du(m->nu);
    dx.setZero(); du.setZero();
    // Predictions
    eigVm xNewNominal(cp.n), xNewPerturbed(cp.n), xNewHW(cp.n), xNewW(cp.n), xNewP(cp.n);
    // Computation time and error vectors
    eigVd tHW(nSample), tW(nSample), tP(nSample), eHW(nSample), eW(nSample), eP(nSample), fCon(nSample);
    tHW.setZero(); tW.setZero(); tP.setZero(); eHW.setZero(); eW.setZero(); eP.setZero(); fCon.setZero();
    eigVd acc_err(nSample); acc_err.setZero();

    // Set random seed
    // std::srand(std::time(0));
    std::srand(0);

    // Test loop
    for( int k=0; k<nSample; k++ )
    {
        std::cout << "\n====================== Sample no: " << k+1 << " ======================\n";
        // Generate random configuration
        if( k>0 )
        {
            qRand = Eigen::VectorXd::Random(m->nq)*5e-1;
            vRand = Eigen::VectorXd::Random(m->nv)*5e-2;
        }
        // Create MuJoCo data
        mjData* d = mj_makeData(m);
        // Set the joint positions to the preset configuration in the model
        if(k==0)
            mju_copy(d->qpos, m->key_qpos, m->nq);
        
        // Copy random variables to MuJoCo data
        mju_copy(d->qpos, qRand.data(), m->nq);
        mju_copy(d->qvel, vRand.data(), m->nv);

        // Run forward dynamics to calculate the bias
        mj_forward(m, d);        
        mju_copy(d->ctrl, d->qfrc_bias, m->nv);

        // Run forward dynamics for full computation
        mj_forward(m, d);
        
        // get initial state & control
        x = cc.getState(d);
        mju_copy(u.data(), d->ctrl, m->nu);

        // Copy MuJoCo data to Pinocchio variables
        mju_copy(q.data(), d->qpos, m->nq);
        mju_copy(v.data(), d->qvel, m->nv);
        tau = u;

        // Print initial states and controls
        if( printMInit )
        {
            printConfig(m, d);
        }
        if( printPInit )
        {
            std::cout<< "Pinocchio state and control before perturbation:\n";
            std::cout << "  q   = " << q.transpose() << "\n";
            std::cout << "  v   = " << v.transpose() << "\n";
            std::cout << "  tau = " << tau.transpose() << "\n\n";
        }

        // Get contact forces projected onto the joint space
        mju_copy(qcon.data(), d->qfrc_constraint, m->nv);
        fCon[k] = qcon.norm();
        
        // Print contact info
        if(printConst)
            printContactInfo(m, d);
        // Get contact forces from MuJoCo and set fext
        PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext;
        fext = getExternalForce(m, d, model);

        // Compare accelerations
        Eigen::VectorXd mj_a(m->nv), pin_a(model.nv);
        mju_copy(mj_a.data(), d->qacc, m->nv);
        pinocchio::aba(model, data, q, v, tau, fext);
        pin_a = data.ddq;
        pinocchio::aba(model, data, q, v, tau+qcon);
        std::cout << "Accelerations:\n" <<
                     "\tMuJoCo:             " << mj_a.transpose() <<
                     "\n\tPinocchio w/o fext: " << data.ddq.transpose() <<
                     "\n\tPinocchio w/ fext:  " << pin_a.transpose() << "\n\n";
        acc_err(k) = (mj_a - pin_a).norm();

        // Calculate derivatives with MuJoCo hardWorker
        auto tHWStart = std::chrono::system_clock::now();
        nd.linDyn(d, u, FxHW.data(), FuHW.data(), compensateBias);
        auto tHWEnd = std::chrono::system_clock::now();
        tHW(k) = std::chrono::duration<double>(tHWEnd-tHWStart).count();
        if( printTime )
            std::cout << "\nINFO: MuJoCo hardWorker took " << tHW(k) << " s\n\n";

        // Calculate derivatives with MuJoCo worker
        auto tMjStart = std::chrono::system_clock::now();
        nd.worker(d, deriv);
        auto tMjEnd = std::chrono::system_clock::now();
        tW(k) = std::chrono::duration<double>(tMjEnd-tMjStart).count();
        if( printTime )
            std::cout << "\nINFO: MuJoCo worker took " << tW(k) << " s\n\n";
        // get derivatives of acceleration
        mju_copy(da_dq.data(), deriv, m->nv*m->nv);
        mju_copy(da_dv.data(), deriv+m->nv*m->nv, m->nv*m->nv);
        mju_copy(da_df.data(), deriv+2*m->nv*m->nv, m->nv*m->nv);
        // build derivative matrices
        FxW.setZero();
        FxW.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*cp.tc*da_dq;
        FxW.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*cp.tc*da_dv;
        FxW.bottomLeftCorner(m->nv, m->nv)  = cp.tc*da_dq;
        FxW.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*da_dv;
        FuW.setZero();
        FuW.topRows(m->nv)    = cp.tc*cp.tc*da_df;
        FuW.bottomRows(m->nv) = cp.tc*da_df;

        // Calculate derivatives with Pinocchio
        auto tPinStart = std::chrono::system_clock::now();
        if(fext_flag==0)
            pinocchio::computeABADerivatives(model, data, q, v, tau);
        else if(fext_flag==1)
            pinocchio::computeABADerivatives(model, data, q, v, tau, fext);
        else if(fext_flag==2)
            pinocchio::computeABADerivatives(model, data, q, v, tau+qcon);
        auto tPinEnd = std::chrono::system_clock::now();
        tP(k) = std::chrono::duration<double>(tPinEnd-tPinStart).count();
        if( printTime )
            std::cout << "\nINFO: Pinocchio took " << tP(k) << " s\n\n";
        // build derivative matrices
        FxP.setZero();
        FxP.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*cp.tc*data.ddq_dq;
        FxP.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*cp.tc*data.ddq_dv;
        FxP.bottomLeftCorner(m->nv, m->nv)  = cp.tc*data.ddq_dq;
        FxP.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*data.ddq_dv;
        FuP.setZero();
        FuP.topRows(m->nv)    = cp.tc*cp.tc*data.Minv;
        FuP.bottomRows(m->nv) = cp.tc*data.Minv;

        // Print derivative matrices
        if( printDeriv )
        {
            std::cout << "FxHW:\n" << FxHW << "\n\n";
            std::cout << "FxW:\n" << FxW << "\n\n";
            std::cout << "FxP:\n" << FxP << "\n\n";
            std::cout << "FuHW:\n" << FuHW << "\n\n";
            std::cout << "FuW:\n" << FuW << "\n\n";
            std::cout << "FuP:\n" << FuP << "\n\n";
        }

        // Nominal trajectory
        mjData* dNominal = mj_makeData(m);
        copyData(m, d, dNominal);
        mj_forward(m, dNominal);
        if( printTraj )
        {
            std::cout << "Nominal trajectory:";
            printConfig(m, dNominal);
        }
        cc.takeStep(dNominal, u, false, compensateBias);
        if( printTraj )
        { printConfig(m, dNominal); }
        xNewNominal = cc.getState(dNominal);
        mj_deleteData(dNominal);

        // Perturbation
        if(pos_pert)
            dx.head(m->nv) = Eigen::VectorXd::Random(m->nv)*1e-2;
        if(vel_pert)
            dx.tail(m->nv) = Eigen::VectorXd::Random(m->nv)*1e-2;
        if(tau_pert)
            du = Eigen::VectorXd::Random(m->nu)*1e-2;
        // print perturbation
        if( printPert )
        {
            std::cout << "Perturbation:\n";
            std::cout << "  x: " << x.transpose() << "\n  dq: " << dx.head(m->nv).transpose() << "\n";
            std::cout << "  dv: " << dx.tail(m->nv).transpose() << "\n";
            std::cout << "  xnew: " << (x+dx).transpose() << "\n\n";
            std::cout << "  u: " << u.transpose() << "\n  du: " << du.transpose() << "\n";
            std::cout << "  unew: " << (u+du).transpose() << "\n";
        }
        // apply perturbation
        x += dx;
        u += du;

        // Perturbed trajectory
        mjData* dPerturbed = mj_makeData(m);
        copyData(m, d, dPerturbed);
        mju_copy(dPerturbed->ctrl, u.data(), m->nu);
        // mju_copy(dPerturbed->qfrc_applied, u.data(), m->nv);
        mju_copy(dPerturbed->qpos, x.head(m->nv).data(), m->nv);
        mju_copy(dPerturbed->qvel, x.tail(m->nv).data(), m->nv);
        mj_forward(m, dPerturbed);
        // Print contact info
        if(printConst)
            printContactInfo(m, d);
        // Integrate w/ perturbations
        cc.takeStep(dPerturbed, u, false, compensateBias);
        // Get the perturbed state and delete data
        xNewPerturbed = cc.getState(dPerturbed);
        mj_deleteData(dPerturbed);

        // Predictions
        xNewHW.setZero();   xNewW.setZero();   xNewP.setZero();
        xNewHW = xNewNominal + FxHW*dx + FuHW*du;
        xNewW  = xNewNominal + FxW*dx + FuW*du;
        xNewP  = xNewNominal + FxP*dx + FuP*du;
        eHW(k) = (xNewPerturbed-xNewHW).norm();
        eW(k)  = (xNewPerturbed-xNewW).norm();
        eP(k)  = (xNewPerturbed-xNewP).norm();
        if( printPred )
        {
            std::cout << "Perturbation:\n";
            std::cout << "  dq: " << dx.head(m->nv).transpose() << "\n";
            std::cout << "  dv: " << dx.tail(m->nv).transpose() << "\n";
            std::cout << "  du: " << du.transpose() << "\n";
            std::cout << "Nominal next state:\n";
            std::cout << "  pos: " << xNewNominal.head(m->nv).transpose() << "\n";
            std::cout << "  vel: " << xNewNominal.tail(m->nv).transpose() << "\n";
            std::cout << "Perturbed next state:\n";
            std::cout << "  pos: " << xNewPerturbed.head(m->nv).transpose() << "\n";
            std::cout << "  vel: " << xNewPerturbed.tail(m->nv).transpose() << "\n";
            // hardWorker
            std::cout << "hardWorker prediction:\n";
            std::cout << "  pos: " << xNewHW.head(m->nv).transpose() << "\n";
            std::cout << "  vel: " << xNewHW.tail(m->nv).transpose() << "\n";
            std::cout << "  norm(error): " << eHW(k) << "\n";
            std::cout << "  comp. time: " << tHW(k) << " s\n";
            // worker
            std::cout << "worker prediction:\n";
            std::cout << "  pos: " << xNewW.head(m->nv).transpose() << "\n";
            std::cout << "  vel: " << xNewW.tail(m->nv).transpose() << "\n";
            std::cout << "  norm(error): " << eW(k) << "\n";
            std::cout << "  comp. time: " << tW(k) << " s\n";
            // Pinocchio
            std::cout << "pinocchio prediction:\n";
            std::cout << "  pos: " << xNewP.head(m->nv).transpose() << "\n";
            std::cout << "  vel: " << xNewP.tail(m->nv).transpose() << "\n";
            std::cout << "  norm(error): " << eP(k) << "\n";
            std::cout << "  comp. time: " << tP(k) << " s\n";
        }
        // Delete data
        mj_deleteData(d);
    }
    std::cout << "\n================================================================================\n\n";
    if(fext_flag==0)
        std::cout << "External forces were not taken into account in Pinocchio derivatives.\n";
    else if(fext_flag==1)
        std::cout << "External forces were taken into account in Pinocchio derivatives.\n";
    else if(fext_flag==2)
        std::cout << "Joint-space contact forces were added to Pinocchio torques.\n";
    std::cout << "INFO: Test done.\n\nTest summary:\n";
    for( int k=0; k<nSample; k++ )
    {
        if( k%10 == 0 )
        {
            printf("%-10s%-10s%-14s%-36s%-36s\n",
                   "Sample","Contact","Accel.","Prediction Error","Computation Time [s]");
            printf("%-10s%-10s%-14s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "#","force","error","hardWorker","worker","Pinocchio","hardWorker","worker","Pinocchio");
        }
        printf("%-10d%-10.6g%-14.6g%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g\n",
               k+1,fCon(k),acc_err(k),eHW(k),eW(k),eP(k),tHW(k),tW(k),tP(k));
    }

    printf("\nStatistics for %d samples:\n",nSample);
    printf("%-16s%-16s%-16s%-16s\n",
           "Function:","hardWorker","worker","Pinocchio");
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Mean error:",eHW.mean(),eW.mean(),eP.mean());
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Min. error:",eHW.minCoeff(),eW.minCoeff(),eP.minCoeff());
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Max. error:",eHW.maxCoeff(),eW.maxCoeff(),eP.maxCoeff());
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Mean time:",tHW.mean(),tW.mean(),tP.mean());
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Min. time:",tHW.minCoeff(),tW.minCoeff(),tP.minCoeff());
    printf("%-16s%-16.6g%-16.6g%-16.6g\n",
           "Max. time:",tHW.maxCoeff(),tW.maxCoeff(),tP.maxCoeff());


    // Shut down MuJoCo
    mju_free(deriv);
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}