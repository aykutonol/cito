#include <gtest/gtest.h>
#include "cito/scvx.h"

// MuJoCo model and data
mjModel *m = NULL;
mjData *d = NULL;
int mj_activated = 0;

// YAML file for parameters
YAML::Node params;

// Assert model and data are not empty
TEST(CitoTests, test_mujoco_init)
{
    ASSERT_TRUE(mj_activated > 0);
    ASSERT_TRUE(m != NULL);
    ASSERT_TRUE(d != NULL);
}

// Assert params file loaded
TEST(CitoTests, test_params_load)
{
    ASSERT_FALSE(params.IsNull());
}

// Assert skew-symmetric conversion working properly
TEST(CitoTests, test_skew)
{
    // Instantiate a Params object
    Params cp(m);
    // Create a test vector with an expected skew-symmetric representation
    Eigen::Vector3d v(1., 2., 3.);
    Eigen::Matrix3d Vhat_expect;
    Vhat_expect.setZero();
    Vhat_expect(0, 1) = -3.;
    Vhat_expect(0, 2) = 2.;
    Vhat_expect(1, 2) = -1.;
    Vhat_expect(1, 0) = 3.;
    Vhat_expect(2, 0) = -2.;
    Vhat_expect(2, 1) = 1.;
    // Apply the util function for skew-symmetric conversion
    Eigen::Matrix3d Vhat = cp.skew(v);
    // Compare the resulting and expected matrices
    ASSERT_TRUE(Vhat.isApprox(Vhat_expect));
}

// Assert contact normal Jacobian calculation working properly
TEST(CitoTests, test_evalNormalJac)
{
    // Instantiate a Params object
    Params cp(m);
    // Set the joint configuration and evaluate the forward dynamics
    mju_copy(d->qpos, m->key_qpos, m->nq);
    srand((unsigned int)time(0));
    Eigen::Vector4d quat_rand = Eigen::MatrixXd::Random(4, 1);
    // quat_rand << 0.842, 0.408, 0.342, 0.092;
    // quat_rand << .966, .0, .0, .259;
    quat_rand = quat_rand / quat_rand.norm();
    mju_copy4(d->qpos + 3, quat_rand.data());
    std::cout << "qrand: " << quat_rand.transpose() << ", norm: " << quat_rand.norm() << "\nEuler: " << cp.quat2Euler(quat_rand).transpose() << "\n\n";
    mj_forward(m, d);
    // Evaluate the Jacobian numerically
    Eigen::Vector3d n, np;
    Eigen::MatrixXd dn_dw(3, m->nv);
    double eps = 1e-6;
    for (int pair = 0; pair < cp.nPair; pair++)
    {
        mju_rotVecMat(n.data(), cp.unit_x, d->site_xmat + 9 * cp.sites[pair][1]);
        for (int i = 0; i < m->nv; i++)
        {
            mjData *dCopy = mj_makeData(m);
            mj_copyData(dCopy, m, d);
            int jID = m->dof_jntid[i];
            if (cp.quatAdr[i] >= 0)
            {
                mjtNum angvel[3] = {0, 0, 0};
                angvel[cp.dofAdr[i]] = eps;
                mju_quatIntegrate(dCopy->qpos + cp.quatAdr[i], angvel, 1);
            }
            else
                dCopy->qpos[m->jnt_qposadr[jID] + i - m->jnt_dofadr[jID]] += eps;
            // get the perturbed state
            mj_forward(m, dCopy);
            mju_rotVecMat(np.data(), cp.unit_x, dCopy->site_xmat + 9 * cp.sites[pair][1]);
            // undo perturbation
            mju_copy(dCopy->qpos, d->qpos, m->nq);
            // compute column i of dn/dq
            dn_dw.col(i) = (np - n) / eps;
            mj_deleteData(dCopy);
        }
        // Evaluate the analytic Jacobian
        Eigen::Vector4d quat_site;
        mju_mat2Quat(quat_site.data(), d->site_xmat + 9 * cp.sites[pair][1]);
        Eigen::Matrix3d dn_dw_a = cp.evalNormalJac(quat_site, pair);
        std::cout << "pair: " << pair << ", n: " << n.transpose() << ", Euler: " << cp.quat2Euler(quat_site).transpose() << "\n  numeric dn_dw:\n"
                  << dn_dw.block(0, 3, 3, 3) << "\n  analytic dn_dw:\n"
                  << dn_dw_a << "\n\n";
        // Compare the resulting and expected matrices
        ASSERT_TRUE((dn_dw.block(0, 3, 3, 3) - dn_dw_a).cwiseAbs().maxCoeff() < 1e-3);
    }
}

int main(int argc, char **argv)
{
    // Initialize testing
    testing::InitGoogleTest(&argc, argv);

    // Model file
    params = YAML::LoadFile(paths::workspaceDir + "/src/cito/config/params.yaml");
    std::string modelPathStr = paths::workspaceDir + "/src/cito/model/" + params["model"].as<std::string>();
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";

    // Activate MuJoCo
    const char *mjKeyPath = std::getenv("MJ_KEY");
    mj_activated = mj_activate(mjKeyPath);
    // Load the model
    if (strlen(modelPath) > 4 && !strcmp(modelPath + strlen(modelPath) - 4, ".mjb"))
    {
        m = mj_loadModel(modelPath, NULL);
    }
    else
    {
        m = mj_loadXML(modelPath, NULL, NULL, 0);
    }
    if (!m)
    {
        mju_error("Cannot load the model");
    }
    // Create data
    d = mj_makeData(m);

    // Run and return tests
    return RUN_ALL_TESTS();
}