#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

int ndof=6, pos_off=0, vel_off=0;

// Compensate bias term
double compensateBias = 1.0;

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

/// Flags
bool printDeriv = false;
bool printMInit = false;
bool printPInit = false;
bool printPred  = false;
bool printPert  = false;
bool printTraj  = false;
bool printTime  = false;

int main(int argc, char const *argv[]) {
    /// Parse command line arguments
    int nSample = 50, fext_flag=2;
    if(argc>1)
    {
        nSample = atoi(argv[1]);
        if(argc>2)
        {
            fext_flag = atoi(argv[2]);
        }
    }
    /// Initialize MuJoCo
    // Activate MuJoCo
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
    /// Initialize Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/ur3e.urdf", model);
    pinocchio::Data data(model);
    /// Initialize CITO objects
    CitoParams  cp(m);
    CitoNumDiff nd(m);
    CitoControl cc(m);
    /// Initialize variables
    // state and control vectors for MuJoCo
    eigVm x, u;
    x.resize(cp.n); u.resize(cp.m);
    x.setZero();    u.setZero();
    // state and control matrices for Pinocchio
    eigVd q, v, tau, qcon;
    q.resize(model.nq); v.resize(model.nv); tau.resize(model.nv); qcon.resize(model.nv);
    q.setZero();        v.setZero();        tau.setZero();        qcon.setZero();
    PINOCCHIO_ALIGNED_STD_VECTOR(pinocchio::Force) fext((size_t)model.njoints, pinocchio::Force::Zero());
    // derivative matrices
    eigMd da_dq, da_dv, da_du;
    da_dq.resize(m->nv, m->nv); da_dv.resize(m->nv, m->nv); da_du.resize(m->nv, m->nu);
    eigMd FxHW, FuHW, FxW, FuW, FxP, FuP;
    FxHW.resize(cp.n, cp.n);    FxW.resize(cp.n, cp.n);     FxP.resize(cp.n, cp.n);
    FuHW.resize(cp.n, cp.m);    FuW.resize(cp.n, cp.m);     FuP.resize(cp.n, cp.m);
    // random configurations
    eigVd qRand, vRand, uRand;
    // perturbations
    eigVm dx(cp.n), du(m->nu);
    dx.setZero();       du.setZero();
    // predictions
    eigVm xNewNominal(cp.n), xNewPerturbed(cp.n), xNewHW(cp.n), xNewW(cp.n), xNewP(cp.n);
    // computation time and error vectors
    eigVd tHW(nSample), tW(nSample), tP(nSample), eHW(nSample), eW(nSample), eP(nSample), fCon(nSample);
    tHW.setZero(); tW.setZero(); tP.setZero(); eHW.setZero(); eW.setZero(); eP.setZero(); fCon.setZero();

    /// Set random seed
    // std::srand(std::time(0));
    std::srand(0);

    /// Test loop
    for( int k=0; k<nSample; k++ )
    {
        std::cout << "\n====================== Sample no: " << k+1 << " ======================\n";
        /// Generate random configuration
        if( k==0 )
        {
            // qRand = Eigen::VectorXd::Zero(m->nq);
            vRand = Eigen::VectorXd::Zero(m->nv);
            // uRand = Eigen::VectorXd::Zero(m->nu);
        }
        else
        {
            // qRand = Eigen::VectorXd::Random(m->nq)*2;
            vRand = Eigen::VectorXd::Random(m->nv)*1;
            // uRand = Eigen::VectorXd::Random(m->nu)*1;
        }
        // Create MuJoCo data
        mjData* d = mj_makeData(m);
        mjData* dTemp = mj_makeData(m);
        // Set the joint positions to the preset configuration in the model
        mju_copy(d->qpos, m->key_qpos, m->nq);
        mj_forward(m, d);
        // Copy random variables to MuJoCo data
        // mju_copy(d->qpos, qRand.data(), m->nq);
        mju_copy(d->qvel, vRand.data(), m->nv);
        // mju_copy(d->ctrl, uRand.data(), m->nu);
        // mju_copy(d->ctrl, d->qfrc_bias, m->nu);
        mju_copy(d->qfrc_applied, d->qfrc_bias, m->nv);
        
        // Take warm-up steps to ensure making contacts
        int n_warmup = 5;
        for(int i=0; i<n_warmup; i++)
        {
            mj_step(m, d);
        }
        // get initial state & control
        x = cc.getState(d);
        mju_copy(u.data(), d->qfrc_applied, m->nv);

        if( printMInit )
        {
            printConfig(m, d);
        }

        // Copy MuJoCo data to Pinocchio variables
        mju_copy(q.data(), d->qpos+pos_off, model.nq);
        mju_copy(v.data(), d->qvel+vel_off, model.nv);
        mju_copy(tau.data(), d->qfrc_applied+vel_off, model.nv);
        mju_copy(qcon.data(), d->qfrc_constraint, m->nv);
        fCon[k] = qcon.norm();

        // Constraints
        std::cout << "\n\nConstraints:\n";
        for(int i=0; i<d->nefc; i++)
        {
            std::cout << "\tConstraint " << i << ", type: " << d->efc_type[i] <<
                         ", force: " << d->efc_force[i] <<
                         ", state: " << d->efc_state[i] <<
                         ", id: " << d->efc_id[i] << "\n";
        }
        // Contacts
        eigVd hcon(6); hcon.setZero();
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

        // Set external force in joint space
        if(fext_flag>0)
        {
            pinocchio::Force::Vector3 tau_contact = pinocchio::Force::Vector3::Zero();
            for(int i=1; i<model.njoints; i++)
            {
                // std::cout << "MuJoCo joint axis: " << m->jnt_axis[(i-1)*3] << m->jnt_axis[(i-1)*3+1] << m->jnt_axis[(i-1)*3+2] << "\n";
                // std::cout << model.joints[i] << "\n";
                for(int j=0; j<3; j++)
                {
                    if(fext_flag==1)
                        tau_contact[j] =  m->jnt_axis[(i-1)*3+j]*d->qfrc_constraint[i-1];
                    else if(fext_flag==2)
                        tau_contact[j] =  -m->jnt_axis[(i-1)*3+j]*d->qfrc_constraint[i-1];
                }
                fext[i].angular(tau_contact);
                // std::cout << "tau_contact: " << tau_contact.transpose() << "\n";
                // std::cout << fext[i] << "\n";
            }
        }

        if( printPInit )
        {
            std::cout<< "Pinocchio state and control before perturbation:\n";
            std::cout << "  q   = " << q.transpose() << "\n";
            std::cout << "  v   = " << v.transpose() << "\n";
            std::cout << "  tau = " << tau.transpose() << "\n\n";
        }

        /// Calculate derivatives with MuJoCo hardWorker
        auto tHWStart = std::chrono::system_clock::now();
        nd.linDyn(d, u, FxHW.data(), FuHW.data(), compensateBias);
        auto tHWEnd = std::chrono::system_clock::now();
        tHW(k) = std::chrono::duration<double>(tHWEnd-tHWStart).count();
        if( printTime )
            std::cout << "\nINFO: MuJoCo hardWorker took " << tHW(k) << " s\n\n";

        /// Calculate derivatives with MuJoCo worker
        auto tMjStart = std::chrono::system_clock::now();
        nd.worker(d, deriv);
        auto tMjEnd = std::chrono::system_clock::now();
        tW(k) = std::chrono::duration<double>(tMjEnd-tMjStart).count();
        if( printTime )
            std::cout << "\nINFO: MuJoCo worker took " << tW(k) << " s\n\n";
        // get derivatives of acceleration
        mju_copy(da_dq.data(), deriv, m->nv*m->nv);
        mju_copy(da_dv.data(), deriv+m->nv*m->nv, m->nv*m->nv);
        mju_copy(da_du.data(), deriv+2*m->nv*m->nv, m->nv*m->nu);
        // build derivative matrices
        FxW.setZero();
        FxW.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxW.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxW.bottomLeftCorner(m->nv, m->nv)  = cp.tc*da_dq;
        FxW.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*da_dv;
        FuW.setZero();
        FuW.bottomRows(m->nv) = cp.tc*da_du;

        /// Calculate derivatives with Pinocchio
        auto tPinStart = std::chrono::system_clock::now();
        if(fext_flag==0)
            pinocchio::computeABADerivatives(model, data, q, v, tau);
        else if(fext_flag==1 || fext_flag==2)
            pinocchio::computeABADerivatives(model, data, q, v, tau, fext);
        else if(fext_flag==3)
            pinocchio::computeABADerivatives(model, data, q, v, tau+qcon);
        auto tPinEnd = std::chrono::system_clock::now();
        tP(k) = std::chrono::duration<double>(tPinEnd-tPinStart).count();
        if( printTime )
            std::cout << "\nINFO: Pinocchio took " << tP(k) << " s\n\n";
        // build derivative matrices
        FxP.setZero();
        FxP.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxP.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxP.bottomLeftCorner(m->nv, m->nv)  = cp.tc*data.ddq_dq;
        FxP.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*data.ddq_dv;
        FuP.setZero();
        FuP.bottomRows(m->nv) = cp.tc*data.Minv;

        /// Print derivative matrices
        if( printDeriv )
        {
            std::cout << "FxHW:\n" << FxHW << "\n\n";
            std::cout << "FxW:\n" << FxW << "\n\n";
            std::cout << "FxP:\n" << FxP << "\n\n";
            std::cout << "FuHW:\n" << FuHW << "\n\n";
            std::cout << "FuW:\n" << FuW << "\n\n";
            std::cout << "FuP:\n" << FuP << "\n\n";
        }

        /// Nominal trajectory
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

        /// Perturbation
        // generate random perturbation
        if(k==0)
        {
            dx = Eigen::VectorXd::Zero(cp.n);
            du = Eigen::VectorXd::Zero(cp.m);
        }
        else
        {
            dx = Eigen::VectorXd::Random(cp.n)*1e-2;
            du = Eigen::VectorXd::Random(cp.m)*1e-2;
        }
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

        /// Perturbed trajectory
        mjData* dPerturbed = mj_makeData(m);
        copyData(m, d, dPerturbed);
        // mju_copy(dPerturbed->ctrl, u.data(), m->nu);
        mju_copy(dPerturbed->qfrc_applied, u.data(), m->nv);
        mju_copy(dPerturbed->qpos, x.head(m->nv).data(), m->nv);
        mju_copy(dPerturbed->qvel, x.tail(m->nv).data(), m->nv);
        mj_forward(m, dPerturbed);
        if( printTraj )
        {
            std::cout << "Perturbed trajectory:";
            printConfig(m, dPerturbed);
        }
        cc.takeStep(dPerturbed, u, false, compensateBias);
        if( printTraj )
        { printConfig(m, dPerturbed); }
        xNewPerturbed = cc.getState(dPerturbed);
        mj_deleteData(dPerturbed);

        /// Predictions
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
            std::cout << "Actual next state:\n";
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
        mj_deleteData(dTemp);
    }
    std::cout << "\n================================================================================\n\n";
    std::cout << "INFO: Test done.\n\nTest summary:\n";
    for( int k=0; k<nSample; k++ )
    {
        if( k%10 == 0 )
        {
            printf("%-12s%-12s%-36s%-36s\n",
                   "Sample","Contact","Prediction Error","Computation Time [s]");
            printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "#","force","hardWorker","worker","Pinocchio","hardWorker","worker","Pinocchio");
        }
        printf("%-12d%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g\n",
               k+1,fCon(k),eHW(k),eW(k),eP(k),tHW(k),tW(k),tP(k));
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