#include <chrono>

#include "cito_scvx.h"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"

mjtNum* deriv = 0;          // dynamics derivatives (6*nv*nv):
int nwarmup = 3;            // center point repetitions to improve warmstart
double eps = 1e-6;          // finite-difference epsilon

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
    std::cout << "qacc: ";
    mju_printMat(d->qacc, 1, m->nv);
    std::cout << "tau: ";
    mju_printMat(d->ctrl, 1, m->nu);
}

void copyData(mjModel* m, mjData* dmain, mjData* d)
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
    mjData* dtemp = mj_makeData(m);

    d->qvel[0] = 1;

    mj_forward(m, d);
    showConfig(m, d);
    // Allocate derivatives
    deriv = (mjtNum*) mju_malloc(3*sizeof(mjtNum)*m->nv*m->nv);
    // Calculate derivatives
    auto tMjStart = std::chrono::system_clock::now();
    worker(m, d, dtemp);
    auto tMjEnd = std::chrono::system_clock::now();
    std::cout << "\nINFO: MuJoCo worker took " << std::chrono::duration<double>(tMjEnd-tMjStart).count() << " s\n\n";

    eigMm da_dq, da_dv, da_df;
    da_dq.resize(m->nv, m->nv); da_dv.resize(m->nv, m->nv); da_df.resize(m->nv, m->nu);
    mju_copy(da_dq.data(), deriv, m->nv*m->nv);
    mju_copy(da_dv.data(), deriv+m->nv*m->nv, m->nv*m->nv);
    mju_copy(da_df.data(), deriv+2*m->nv*m->nv, m->nv*m->nu);

    std::cout << "dqacc/dqpos:\n" << da_dq << "\n----------------\n";
    std::cout << "dqacc/dqvel:\n" << da_dv << "\n----------------\n";
    std::cout << "dqacc/dtau:\n"  << da_df << "\n----------------\n";


    // Pinocchio
    pinocchio::Model model;
    pinocchio::urdf::buildModel("/home/aykut/Development/ur_ws/src/universal_robot/ur_e_description/urdf/ur3e.urdf", model);
    pinocchio::Data data(model);

    eigVd q   = pinocchio::neutral(model);
    eigVd v   = Eigen::VectorXd::Zero(model.nv);
    eigVd a   = Eigen::VectorXd::Zero(model.nv);
    eigVd tau = Eigen::VectorXd::Zero(model.nv);

    for( int i=0; i<model.nv; i++ )
    {
        q[i]   = d->qpos[i];
        v[i]   = d->qvel[i];
        a[i]   = d->qacc[i];
        tau[i] = d->ctrl[i];
    }

    std::cout<< "Pinocchio model name: " << model.name << "\n";
    std::cout << "q   = " << q.transpose() << "\n";
    std::cout << "v   = " << v.transpose() << "\n";
    std::cout << "a   = " << a.transpose() << "\n";
    std::cout << "tau = " << tau.transpose() << "\n";

    auto tPinStart = std::chrono::system_clock::now();
    pinocchio::computeABADerivatives(model, data, q, v, tau);
    auto tPinEnd = std::chrono::system_clock::now();
    std::cout << "\nINFO: Pinocchio took " << std::chrono::duration<double>(tPinEnd-tPinStart).count() << " s\n\n";

    std::cout << "dqacc/dqpos:\n" << data.ddq_dq << "\n\n";
    std::cout << "dqacc/dqvel:\n" << data.ddq_dv << "\n\n";
    std::cout << "dqacc/dtau: \n" << data.Minv << "\n\n";

    // MuJoCo prediction
    CitoParams  cp(m);
    CitoNumDiff nd(m);
    CitoControl cc(m);
    eigMm U;        U = Eigen:: MatrixXd::Zero(cp.m, 1);
    eigMm FxHW, FuHW;
    FxHW.resize(cp.n, cp.n); FuHW.resize(cp.n, cp.m);
    std::cout << "n = " << cp.n << ", m = " << cp.m << ", N: " << cp.N << "\n";
    double tM = cp.dt*cp.ndpc;
    std::cout << "dt = " << cp.dt << ", ndpc = " << cp.ndpc << ", tM = " << tM << "\n";

    // hardWorkder derivatives
    auto tHWStart = std::chrono::system_clock::now();
    nd.linDyn(d, U, FxHW.data(), FuHW.data());
    auto tHWEnd = std::chrono::system_clock::now();
    std::cout << "\nINFO: MuJoCo hardWorker took " << std::chrono::duration<double>(tHWEnd-tHWStart).count() << " s\n\n";


    // worker derivatives
    eigMm FxW; FxW.resize(cp.n, cp.n); FxW.setZero();
    FxW.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxW.topRightCorner(m->nv, m->nv)    = tM*Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxW.bottomLeftCorner(m->nv, m->nv)  = tM*da_dq;
    FxW.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + tM*da_dv;

    eigMm FuW; FuW.resize(cp.n, cp.m); FuW.setZero();
    FuW.bottomRows(m->nv) = tM*da_df;


    // pinocchio derivatives
    eigMm FxP; FxP.resize(cp.n, cp.n); FxP.setZero();
    FxP.topLeftCorner(m->nv, m->nv)     = Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxP.topRightCorner(m->nv, m->nv)    = tM*Eigen::MatrixXd::Identity(m->nv, m->nv);
    FxP.bottomLeftCorner(m->nv, m->nv)  = tM*data.ddq_dq;
    FxP.bottomRightCorner(m->nv, m->nv) = Eigen::MatrixXd::Identity(m->nv, m->nv) + tM*data.ddq_dv;

    eigMm FuP; FuP.resize(cp.n, cp.m); FuP.setZero();
    FuP.bottomRows(m->nv) = tM*data.Minv;

    // show derivative matrices
    std::cout << "FxHW:\n" << FxHW << "\n\n";
    std::cout << "FxW:\n" << FxW << "\n\n";
    std::cout << "FxP:\n" << FxP << "\n\n";
    std::cout << "FuHW:\n" << FuHW << "\n";
    std::cout << "FuW:\n" << FuW << "\n\n";
    std::cout << "FuP:\n" << FuP << "\n\n";

    // Nominal trajectory
    mjData* dNominal = mj_makeData(m);
    copyData(m, d, dNominal);

    std::cout << "Nominal trajectory:";
    showConfig(m, dNominal);
    for( int i=0; i<cp.ndpc; i++ )
        mj_step(m, dNominal);
    showConfig(m, dNominal);
    eigVm xNewN = cc.getState(dNominal);

    // Perturbation
    eigVm x, dx, u, du;
    x.resize(cp.n);  dx.resize(cp.n);   dx.setZero();
    u.resize(m->nu); du.resize(m->nu);  du.setZero();

    YAML::Node perturb = YAML::LoadFile(paths::workspaceDir+"/src/cito/config/perturbation.yaml");
    std::vector<double> dqInput = { perturb["dq"].as<std::vector<double>>() };
    std::vector<double> dvInput = { perturb["dv"].as<std::vector<double>>() };
    std::vector<double> duInput = { perturb["du"].as<std::vector<double>>() };
    dx.head(m->nv) = Eigen::Map<Eigen::VectorXd>(dqInput.data(), dqInput.size());
    dx.tail(m->nv) = Eigen::Map<Eigen::VectorXd>(dvInput.data(), dvInput.size());
    du = Eigen::Map<Eigen::VectorXd>(duInput.data(), duInput.size());

    mju_copy(u.data(), d->ctrl, m->nu);
    std::cout << "  u: " << u.transpose() << "\n  du: " << du.transpose() <<  "\n";
    u += du;
    std::cout << "  unew: " << u.transpose() << "\n";

    x = cc.getState(d);
    std::cout << "  x: " << x.transpose() << "\n  dq: " << dx.head(m->nv).transpose() << "\n";
    std::cout << "  dv: " << dx.tail(m->nv).transpose() << "\n";
    x += dx;
    std::cout << "  xnew: " << x.transpose() << "\n\n";


    // Perturbed trajectory
    mjData* dPerturbed = mj_makeData(m);
    copyData(m, d, dPerturbed);
    mju_copy(dPerturbed->ctrl, u.data(), m->nu);
    mju_copy(dPerturbed->qpos, x.head(m->nv).data(), m->nv);
    mju_copy(dPerturbed->qvel, x.tail(m->nv).data(), m->nv);
    // run full computation at center point (usually faster than copying dmain)
    mj_forward(m, dPerturbed);

    std::cout << "Perturbed trajectory:";
    showConfig(m, dPerturbed);
    for( int i=0; i<cp.ndpc; i++ )
    {
        mj_step(m, dPerturbed);
    }
    showConfig(m, dPerturbed);
    eigVm xNewPerturbed = cc.getState(dPerturbed);

    std::cout << "Perturbation:\n";
    std::cout << "  dx: " << dx.transpose() << "\n";
    std::cout << "  du: " << du.transpose() << "\n\n";

    std::cout << "actual next state:\n";
    std::cout << "  pos: " << xNewPerturbed.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewPerturbed.tail(m->nv).transpose() << "\n\n";

    eigVm xNewHW(cp.n); xNewHW.setZero();
    xNewHW = xNewN + FxHW*dx + FuHW*du;
    std::cout << "hardWorker prediction:\n";
    std::cout << "  pos: " << xNewHW.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewHW.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewHW).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tHWEnd-tHWStart).count() << " s\n\n";

    eigVm xNewW(cp.n); xNewW.setZero();
    xNewW = xNewN + FxW*dx + FuW*du;
    std::cout << "worker prediction:\n";
    std::cout << "  pos: " << xNewW.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewW.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewW).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tMjEnd-tMjStart).count() << " s\n\n";

    eigVm xNewP(cp.n); xNewP.setZero();
    xNewP = xNewN + FxP*dx + FuP*du;
    std::cout << "pinocchio prediction:\n";
    std::cout << "  pos: " << xNewP.head(m->nv).transpose() << "\n";
    std::cout << "  vel: " << xNewP.tail(m->nv).transpose() << "\n";
    std::cout << "  norm(error): " << (xNewPerturbed-xNewP).norm() << "\n";
    std::cout << "  comp. time: " << std::chrono::duration<double>(tPinEnd-tPinStart).count() << " s\n\n";

    // Shut down MuJoCo
    mju_free(deriv);
    mj_deleteData(d);
    mj_deleteData(dtemp);
    mj_deleteModel(m);
    mj_deactivate();
    return 0;
}