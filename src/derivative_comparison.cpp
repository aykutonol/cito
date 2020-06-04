#include <chrono>

#include "cito_numdiff.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

void showConfig(mjModel* m, mjData* d)
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
}

void copyData(const mjModel* m, const mjData* dmain, mjData* d)
{
    // Copy state and control from dmain to thread-specific d
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
bool showDeriv = false;
bool showMInit = false;
bool showPInit = true;
bool showPred  = true;
bool showPert  = false;
bool showTraj  = false;
bool showTime  = false;
bool readYAML  = false;

// Compensate bias term
double compensateBias = 0.0;

// Number of samples
int nSample = 50;

int main(int argc, char const *argv[]) {
    // Get the number of samples from the command, if set
    if(argc>1)
    {
        nSample = atoi(argv[1]);
    }
    // Get workspace directory path
    const std::string workspaceDir = std::getenv("CITO_WS");

    // Initialize MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    /// Load model
    std::string mjModelPathStr = paths::workspaceDir + "/src/cito/model/ur3e.xml";
    const char *mjModelPath = mjModelPathStr.c_str();
    mjModel* m = mj_loadXML(mjModelPath, NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    /// Create data
    mjData* d = mj_makeData(m);
    mjData* dTemp = mj_makeData(m);
    /// Allocate memory for derivatives
    mjtNum* deriv = 0;
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);

    /// Initialize CITO objects
    CitoParams  cp(m);
    CitoNumDiff nd(m);
    CitoControl cc(m);
    
    // Initialize Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel(paths::workspaceDir+"/src/cito/model/ur3e.urdf", model);
    pinocchio::Data data(model);

    // Initialize derivative calculation variables
    /// State and control vectors for MuJoCo
    eigVm x, u;
    x.resize(cp.n); u.resize(cp.m);
    x.setZero();    u.setZero();
    /// State and control vectors for Pinocchio
    eigVd q, v, tau;
    q.resize(model.nq); v.resize(model.nv); tau.resize(model.nv);
    q.setZero();        v.setZero();        tau.setZero();
    /// Derivative matrices
    eigMd da_dq, da_dv, da_df;
    da_dq.resize(m->nv, m->nv); da_dv.resize(m->nv, m->nv); da_df.resize(m->nv, m->nu);
    eigMd FxHW, FuHW, FxW, FuW, FxP, FuP;
    FxHW.resize(cp.n, cp.n);    FxW.resize(cp.n, cp.n);     FxP.resize(cp.n, cp.n);
    FuHW.resize(cp.n, cp.m);    FuW.resize(cp.n, cp.m);     FuP.resize(cp.n, cp.m);
    /// Random configuration, velocity, and control vectors
    eigVd qRand, vRand, uRand;
    /// Perturbations
    eigVm dx(cp.n), du(m->nu);
    dx.setZero();       du.setZero();
    /// Predictions
    eigVm xNewNominal(cp.n), xNewPerturbed(cp.n), xNewHW(cp.n), xNewW(cp.n), xNewP(cp.n);
    /// Computation time and error vectors
    eigVd tHW(nSample), tW(nSample), tP(nSample), eHW(nSample), eW(nSample), eP(nSample);
    tHW.setZero(); tW.setZero(); tP.setZero(); eHW.setZero(); eW.setZero(); eP.setZero();
    /// Parse perturbation from file
    if( readYAML )
    {
        YAML::Node perturb = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/perturbation.yaml");
        std::vector<double> dqInput = { perturb["dq"].as<std::vector<double>>() };
        std::vector<double> dvInput = { perturb["dv"].as<std::vector<double>>() };
        std::vector<double> duInput = { perturb["du"].as<std::vector<double>>() };
        dx.head(m->nv) = Eigen::Map<Eigen::VectorXd>(dqInput.data(), dqInput.size());
        dx.tail(m->nv) = Eigen::Map<Eigen::VectorXd>(dvInput.data(), dvInput.size());
        du = Eigen::Map<Eigen::VectorXd>(duInput.data(), duInput.size());
    }

    // Set random seed
    // std::srand(std::time(0));
    std::srand(0);

    // Test loop
    for( int k=0; k<nSample; k++ )
    {
        std::cout << "\n================================================================================\nSample no: " << k+1 << "\n";
        /// Generate random configuration
        if( k==0 )
        {
            qRand = Eigen::VectorXd::Zero(m->nq);
            vRand = Eigen::VectorXd::Zero(m->nv);
            uRand = Eigen::VectorXd::Zero(m->nu);
        }
        else
        {
            qRand = Eigen::VectorXd::Random(m->nq)*2;
            vRand = Eigen::VectorXd::Random(m->nv)*1;
            uRand = Eigen::VectorXd::Random(m->nu)*1;
        }
        // Copy random variables to MuJoCo data
        mju_copy(d->qpos, qRand.data(), m->nq);
        mju_copy(d->qvel, vRand.data(), m->nv);
        mju_copy(d->ctrl, uRand.data(), m->nu);
        // Update control vector
        u = uRand;
        // Evaluate forward dynamics
        mj_forward(m, d);
        // Show random configuration, if opted
        if( showMInit )
        {
            showConfig(m, d);
        }
        // Copy MuJoCo data to Pinocchio variables
        mju_copy(q.data(), d->qpos, m->nq);
        mju_copy(v.data(), d->qvel, m->nv);
        mju_copy(tau.data(), d->ctrl, m->nu);
        if( showPInit )
        {
            std::cout<< "Pinocchio state and control before perturbation:\n";
            std::cout << "  q   = " << q.transpose() << "\n";
            std::cout << "  v   = " << v.transpose() << "\n";
            std::cout << "  tau = " << tau.transpose() << "\n\n";
        }

        // Generate random perturbation
        if( !readYAML )
        {
            dx = Eigen::VectorXd::Random(cp.n)*1e-1;
            du = Eigen::VectorXd::Random(cp.m)*1e-1;
        }

        // Calculate derivatives with MuJoCo hardWorker
        auto tHWStart = std::chrono::system_clock::now();
        nd.linDyn(d, u, FxHW.data(), FuHW.data(), compensateBias);
        auto tHWEnd = std::chrono::system_clock::now();
        tHW(k) = std::chrono::duration<double>(tHWEnd-tHWStart).count();
        if( showTime )
            std::cout << "\nINFO: MuJoCo hardWorker took " << tHW(k) << " s\n\n";

        // Calculate derivatives with MuJoCo worker
        auto tMjStart = std::chrono::system_clock::now();
        nd.worker(d, deriv);
        auto tMjEnd = std::chrono::system_clock::now();
        tW(k) = std::chrono::duration<double>(tMjEnd-tMjStart).count();
        if( showTime )
            std::cout << "\nINFO: MuJoCo worker took " << tW(k) << " s\n\n";
        // Get derivatives of acceleration
        mju_copy(da_dq.data(), deriv, m->nv*m->nv);
        mju_copy(da_dv.data(), deriv+m->nv*m->nv, m->nv*m->nv);
        mju_copy(da_df.data(), deriv+2*m->nv*m->nv, m->nv*m->nu);
        // Build derivative matrices
        FxW.setZero();
        FxW.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxW.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxW.bottomLeftCorner(m->nv, m->nv)  = cp.tc*da_dq;
        FxW.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*da_dv;
        FuW.setZero();
        FuW.bottomRows(m->nv) = cp.tc*da_df;

        // Calculate derivatives with Pinocchio
        auto tPinStart = std::chrono::system_clock::now();
        pinocchio::computeABADerivatives(model, data, q, v, tau);
        auto tPinEnd = std::chrono::system_clock::now();
        tP(k) = std::chrono::duration<double>(tPinEnd-tPinStart).count();
        if( showTime )
            std::cout << "\nINFO: Pinocchio took " << tP(k) << " s\n\n";
        // Build derivative matrices
        FxP.setZero();
        FxP.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxP.topRightCorner(m->nv, m->nv)    = cp.tc*Eigen::MatrixXd::Identity(m->nv, m->nv);
        FxP.bottomLeftCorner(m->nv, m->nv)  = cp.tc*data.ddq_dq;
        FxP.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + cp.tc*data.ddq_dv;
        FuP.setZero();
        FuP.bottomRows(m->nv) = cp.tc*data.Minv;

        // Show derivative matrices
        if( showDeriv )
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
        if( showTraj )
        {
            std::cout << "Nominal trajectory:";
            showConfig(m, dNominal);
        }

        cc.takeStep(dNominal, u, false, compensateBias);

        if( showTraj )
        { showConfig(m, dNominal); }
        xNewNominal = cc.getState(dNominal);
        mj_deleteData(dNominal);

        // Perturbation
        /// Get current state and control
        x = cc.getState(d);
        /// Show perturbation
        if( showPert )
        {
            std::cout << "Perturbation:\n";
            std::cout << "  x: " << x.transpose() << "\n  dq: " << dx.head(m->nv).transpose() << "\n";
            std::cout << "  dv: " << dx.tail(m->nv).transpose() << "\n";
            std::cout << "  xnew: " << (x+dx).transpose() << "\n\n";
            std::cout << "  u: " << u.transpose() << "\n  du: " << du.transpose() << "\n";
            std::cout << "  unew: " << (u+du).transpose() << "\n";
        }
        /// Apply perturbation
        x += dx;
        u += du;

        // Perturbed trajectory
        mjData* dPerturbed = mj_makeData(m);
        copyData(m, d, dPerturbed);
        mju_copy(dPerturbed->ctrl, u.data(), m->nu);
        mju_copy(dPerturbed->qpos, x.head(m->nv).data(), m->nv);
        mju_copy(dPerturbed->qvel, x.tail(m->nv).data(), m->nv);
        mj_forward(m, dPerturbed);
        if( showTraj )
        {
            std::cout << "Perturbed trajectory:";
            showConfig(m, dPerturbed);
        }
        
        cc.takeStep(dPerturbed, u, false, compensateBias);

        if( showTraj )
        { showConfig(m, dPerturbed); }
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
        if( showPred )
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
    }
    std::cout << "\n================================================================================\n\n";
    std::cout << "INFO: Test done.\n\nTest summary:\n";
    for( int k=0; k<nSample; k++ )
    {
        if( k%10 == 0 )
        {
            printf("%-12s%-36s%-36s\n",
                   "Sample","Prediction Error","Computation Time [s]");
            printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n",
                   "#","hardWorker","worker","Pinocchio","hardWorker","worker","Pinocchio");
        }
        printf("%-12d%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g%-12.6g\n",
               k+1,eHW(k),eW(k),eP(k),tHW(k),tW(k),tP(k));
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
    mj_deleteData(d);
    mj_deleteData(dTemp);
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}