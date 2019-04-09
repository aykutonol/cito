#include <chrono>

#include "cito_scvx.h"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

mjtNum* deriv = 0;          // dynamics derivatives (6*nv*nv):
int nwarmup = 3;            // center point repetitions to improve warmstart
double eps = 1e-6;          // finite-difference epsilon

bool showDeriv = true;
bool readYAML  = false;

// worker function for parallel finite-difference computation of derivatives
void worker(const mjModel* m, const mjData* dmain, mjData* d)
{
    int nv = m->nv;

    // allocate stack space for result at center
    mjMARKSTACK
    mjtNum* center = mj_stackAlloc(d, nv);
    mjtNum* warmstart = mj_stackAlloc(d, nv);

    // copy state and control from dmain to d
    d->time = dmain->time;
    mju_copy(d->qpos, dmain->qpos, m->nq);
    mju_copy(d->qvel, dmain->qvel, m->nv);
    mju_copy(d->qacc, dmain->qacc, m->nv);
    mju_copy(d->qacc_warmstart, dmain->qacc_warmstart, m->nv);
    mju_copy(d->qfrc_applied, dmain->qfrc_applied, m->nv);
    mju_copy(d->xfrc_applied, dmain->xfrc_applied, 6*m->nbody);
    mju_copy(d->ctrl, dmain->ctrl, m->nu);

    // run full computation at center point (usually faster than copying dmain)
    mj_forward(m, d);
    // extra solver iterations to improve warmstart (qacc) at center point
    for( int rep=1; rep<nwarmup; rep++ )
        mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

    // select output from forward or inverse dynamics
    mjtNum* output = d->qacc;

    // save output for center point and warmstart (needed in forward only)
    mju_copy(center, output, nv);
    mju_copy(warmstart, d->qacc_warmstart, nv);

    // select target vector and original vector for force or acceleration derivative
    mjtNum* target = d->qfrc_applied;
    const mjtNum* original = dmain->qfrc_applied;

    // finite-difference over force or acceleration: skip = mjSTAGE_VEL
    for( int i=0; i<nv; i++ )
    {
        // perturb selected target
        target[i] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_VEL, 1);

        // undo perturbation
        target[i] = original[i];

        // compute column i of derivative 2
        for( int j=0; j<nv; j++ )
            deriv[2*nv*nv + i + j*nv] = (output[j] - center[j])/eps;
    }

    // finite-difference over velocity: skip = mjSTAGE_POS
    for( int i=0; i<nv; i++ )
    {
        // perturb velocity
        d->qvel[i] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_POS, 1);

        // undo perturbation
        d->qvel[i] = dmain->qvel[i];

        // compute column i of derivative 1
        for( int j=0; j<nv; j++ )
            deriv[nv*nv + i + j*nv] = (output[j] - center[j])/eps;
    }

    // finite-difference over position: skip = mjSTAGE_NONE
    for( int i=0; i<nv; i++ )
    {
        // get joint id for this dof
        int jid = m->dof_jntid[i];

        // get quaternion address and dof position within quaternion (-1: not in quaternion)
        int quatadr = -1, dofpos = 0;
        if( m->jnt_type[jid]==mjJNT_BALL )
        {
            quatadr = m->jnt_qposadr[jid];
            dofpos = i - m->jnt_dofadr[jid];
        }
        else if( m->jnt_type[jid]==mjJNT_FREE && i>=m->jnt_dofadr[jid]+3 )
        {
            quatadr = m->jnt_qposadr[jid] + 3;
            dofpos = i - m->jnt_dofadr[jid] - 3;
        }

        // apply quaternion or simple perturbation
        if( quatadr>=0 )
        {
            mjtNum angvel[3] = {0,0,0};
            angvel[dofpos] = eps;
            mju_quatIntegrate(d->qpos+quatadr, angvel, 1);
        }
        else
            d->qpos[m->jnt_qposadr[jid] + i - m->jnt_dofadr[jid]] += eps;

        // evaluate dynamics, with center warmstart
        mju_copy(d->qacc_warmstart, warmstart, m->nv);
        mj_forwardSkip(m, d, mjSTAGE_NONE, 1);

        // undo perturbation
        mju_copy(d->qpos, dmain->qpos, m->nq);

        // compute column i of derivative 0
        for( int j=0; j<nv; j++ )
            deriv[i + j*nv] = (output[j] - center[j])/eps;
    }

    mjFREESTACK
}

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
    /// Initialize MuJoCo
    // Activate MuJoCo
    const char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // Load model
    mjModel* m;
    m = mj_loadXML("/home/aykut/Development/cito_ws/src/cito/model/ur3e.xml", NULL, NULL, 0);
    if( !m )
        mju_error("Could not load model");
    // Create data
    mjData* d = mj_makeData(m);
    mjData* dTemp = mj_makeData(m);
    // Allocate derivatives
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);
    /// Initialize Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel("/home/aykut/Development/ur_ws/src/universal_robot/ur_e_description/urdf/ur3e.urdf", model);
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
    eigVd q, v, tau;
    q.resize(model.nq); v.resize(model.nv); tau.resize(model.nv);
    q.setZero();        v.setZero();        tau.setZero();
    // time multiplier
    double tM = cp.dt*cp.ndpc;
    // derivative matrices for MuJoCo
    eigMd da_dq, da_dv, da_df;
    da_dq.resize(m->nv, m->nv); da_dv.resize(m->nv, m->nv); da_df.resize(m->nv, m->nu);
    eigMd FxHW, FuHW, FxW, FuW;
    FxHW.resize(cp.n, cp.n);    FuHW.resize(cp.n, cp.m);
    FxW.resize(cp.n, cp.n);     FuW.resize(cp.n, cp.m);
    // derivative matrices for Pinocchio
    eigMd FxP, FuP;
    FxP.resize(cp.n, cp.n);     FuP.resize(cp.n, cp.m);
    /// Read perturbation from YAML
    eigVm dx, du;
    dx.resize(cp.n);    du.resize(m->nu);
    dx.setZero();       du.setZero();
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

    /// Generate random configuration
    eigVd qRand, vRand, uRand;
    qRand = Eigen::VectorXd::Random(m->nq);
    vRand = Eigen::VectorXd::Random(m->nv);
    uRand = Eigen::VectorXd::Random(m->nu);
    // copy random variables to MuJoCo data
    mju_copy(d->qpos, qRand.data(), m->nq);
    mju_copy(d->qvel, vRand.data(), m->nv);
    mju_copy(d->ctrl, uRand.data(), m->nu);
    mj_forward(m, d);
    showConfig(m, d);
    // copy MuJoCo data to Pinocchio variables
    mju_copy(q.data(), d->qpos, m->nq);
    mju_copy(v.data(), d->qvel, m->nv);
    mju_copy(tau.data(), d->ctrl, m->nu);
    std::cout<< "Pinocchio model name: " << model.name << "\n";
    std::cout << "q   = " << q.transpose() << "\n";
    std::cout << "v   = " << v.transpose() << "\n";
    std::cout << "tau = " << tau.transpose() << "\n";

    /// Calculate derivatives with MuJoCo worker
    auto tMjStart = std::chrono::system_clock::now();
    worker(m, d, dTemp);
    auto tMjEnd = std::chrono::system_clock::now();
    auto tW = std::chrono::duration<double>(tMjEnd-tMjStart).count();
    std::cout << "\nINFO: MuJoCo worker took " << tW << " s\n\n";
    // get derivatives of forward dynamics
    mju_copy(da_dq.data(), deriv, m->nv*m->nv);
    mju_copy(da_dv.data(), deriv+m->nv*m->nv, m->nv*m->nv);
    mju_copy(da_df.data(), deriv+2*m->nv*m->nv, m->nv*m->nu);
    // build derivative matrices
    FxW.setZero();
    FxW.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxW.topRightCorner(m->nv, m->nv)    = tM*Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxW.bottomLeftCorner(m->nv, m->nv)  = tM*da_dq;
    FxW.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + tM*da_dv;
    FuW.setZero();
    FuW.bottomRows(m->nv) = tM*da_df;

    /// Calculate derivatives with MuJoCo hardWorker
    auto tHWStart = std::chrono::system_clock::now();
    nd.linDyn(d, u, FxHW.data(), FuHW.data(), 0);
    auto tHWEnd = std::chrono::system_clock::now();
    auto tHW = std::chrono::duration<double>(tHWEnd-tHWStart).count();
    std::cout << "\nINFO: MuJoCo hardWorker took " << tHW << " s\n\n";

    /// Calculate derivatives with Pinocchio
    auto tPinStart = std::chrono::system_clock::now();
    pinocchio::computeABADerivatives(model, data, q, v, tau);
    auto tPinEnd = std::chrono::system_clock::now();
    auto tPin = std::chrono::duration<double>(tPinEnd-tPinStart).count();
    std::cout << "\nINFO: Pinocchio took " << tPin << " s\n\n";
    // build derivative matrices
    FxP.setZero();
    FxP.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxP.topRightCorner(m->nv, m->nv)    = tM*Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxP.bottomLeftCorner(m->nv, m->nv)  = tM*data.ddq_dq;
    FxP.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + tM*data.ddq_dv;
    FuP.setZero();
    FuP.bottomRows(m->nv) = tM*data.Minv;

    /// Show derivative matrices
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
    std::cout << "Nominal trajectory:";
    showConfig(m, dNominal);
    for( int i=0; i<cp.ndpc; i++ )
    {
        mj_step(m, dNominal);
    }
    showConfig(m, dNominal);
    eigVm xNewNominal = cc.getState(dNominal);

    /// Perturbation
    mju_copy(u.data(), d->ctrl, m->nu);
    std::cout << "  u: " << u.transpose() << "\n  du: " << du.transpose() <<  "\n";
    u += du;
    std::cout << "  unew: " << u.transpose() << "\n";

    x = cc.getState(d);
    std::cout << "  x: " << x.transpose() << "\n  dq: " << dx.head(m->nv).transpose() << "\n";
    std::cout << "  dv: " << dx.tail(m->nv).transpose() << "\n";
    x += dx;
    std::cout << "  xnew: " << x.transpose() << "\n\n";


    /// Perturbed trajectory
    mjData* dPerturbed = mj_makeData(m);
    copyData(m, d, dPerturbed);
    mju_copy(dPerturbed->ctrl, u.data(), m->nu);
    mju_copy(dPerturbed->qpos, x.head(m->nv).data(), m->nv);
    mju_copy(dPerturbed->qvel, x.tail(m->nv).data(), m->nv);
    mj_forward(m, dPerturbed);
    std::cout << "Perturbed trajectory:";
    showConfig(m, dPerturbed);
    for( int i=0; i<cp.ndpc; i++ )
    {
        mj_step(m, dPerturbed);
    }
    showConfig(m, dPerturbed);
    eigVm xNewPerturbed = cc.getState(dPerturbed);

    /// Predictions
    std::cout << "Perturbation:\n";
    std::cout << "  dx: " << dx.transpose() << "\n";
    std::cout << "  du: " << du.transpose() << "\n\n";
    std::cout << "Actual next state:\n";
    std::cout << "  pos: " << xNewPerturbed.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewPerturbed.tail(m->nv).transpose() << "\n\n";
    // hardWorker
    eigVm xNewHW(cp.n); xNewHW.setZero();
    xNewHW = xNewNominal + FxHW*dx + FuHW*du;
    std::cout << "hardWorker prediction:\n";
    std::cout << "  pos: " << xNewHW.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewHW.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewHW).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tHWEnd-tHWStart).count() << " s\n\n";
    // worker
    eigVm xNewW(cp.n); xNewW.setZero();
    xNewW = xNewNominal + FxW*dx + FuW*du;
    std::cout << "worker prediction:\n";
    std::cout << "  pos: " << xNewW.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewW.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewW).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tMjEnd-tMjStart).count() << " s\n\n";
    // Pinocchio
    eigVm xNewP(cp.n); xNewP.setZero();
    xNewP = xNewNominal + FxP*dx + FuP*du;
    std::cout << "pinocchio prediction:\n";
    std::cout << "  pos: " << xNewP.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewP.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewP).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tPinEnd-tPinStart).count() << " s\n\n";

    // Shut down MuJoCo
    mju_free(deriv);
    mj_deleteData(d);
    mj_deleteData(dTemp);
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}