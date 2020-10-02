#include "cito/params.h"
#include <fcl/fcl.h>

// Get site pose from MuJoCo data
fcl::Transform3d getSiteTransform(mjData* d, const int& site_id) {
    fcl::Transform3d tf;
    // Get the position
    mju_copy3(tf.translation().data(), d->site_xpos+3*site_id);
    // Get the rotation
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R;
    mju_copy(R.data(), d->site_xmat+9*site_id, 9);
    tf.linear() = R;
    return tf;
}

// Create an FCL collision object given a MuJoCo model, site no, and transform
std::shared_ptr<fcl::CollisionGeometryd> createCollGeom(const mjModel* m, const int& site_id, 
                                                        const fcl::Transform3d& tf) {
    std::shared_ptr<fcl::CollisionGeometryd> geom = NULL;
    if( m->site_type[site_id]==mjGEOM_SPHERE )
    {
        printf("Creating a spherical collision object for site %d on robot.\n", site_id);
        std::cout << "position: " << tf.translation().transpose() << "\n" <<
                     "rotation:\n" << tf.linear() << "\n";
        geom = std::make_shared<fcl::Sphered>(m->site_size[site_id*3]);
    }
    else if( m->site_type[site_id]==mjGEOM_CYLINDER )
    {
        printf("Creating a cylindirical collision object for site %d on robot.\n", site_id);
        std::cout << "position: " << tf.translation().transpose() << "\n" <<
                     "rotation:\n" << tf.linear() << "\n";
        geom = std::make_shared<fcl::Cylinderd>(m->site_size[site_id*3+1]*2,
                                                m->site_size[site_id*3+0]);
    }
    else if( m->site_type[site_id]==mjGEOM_BOX )
    {
        printf("Creating a prismatic collision object for site %d on robot.\n", site_id);
        std::cout << "position: " << tf.translation().transpose() << "\n" <<
                     "rotation:\n" << tf.linear() << "\n";
        geom = std::make_shared<fcl::Boxd>(m->site_size[site_id*3+0]*2, 
                                           m->site_size[site_id*3+1]*2,
                                           m->site_size[site_id*3+2]*2);
    }
    else
    {
        printf("\033[0;31mCannot add site %d b/c geometry type %d is not implemented.\033[0m\n",
               site_id, m->site_type[site_id]);
    }
    return geom;
}

int main(int argc, char* argv[])
{
    // Create distance request & result
    fcl::detail::GJKSolver_libccd<double> GJKSolver;
    fcl::DistanceRequestd distReq;
    fcl::DistanceResultd distRes;
    // Set and print distance request options
    distReq.gjk_solver_type = fcl::GJKSolverType::GST_LIBCCD;
    distReq.enable_nearest_points = true;
    distReq.enable_signed_distance = true;
    std::cout << "\nDistance request options:\n" <<
                 "\ttolerance: " << distReq.distance_tolerance << "\n" <<
                 "\tabs error: " << distReq.abs_err << "\n" <<
                 "\trel error: " << distReq.rel_err << "\n\n";

    // Activate MuJoCo
    char* mjKeyPath = std::getenv("MJ_KEY");
    mj_activate(mjKeyPath);
    // MuJoCo model and data
    mjModel *m = NULL;
    mjData  *d = NULL;
    // Model file
    std::string workspaceDir = std::getenv("CITO_WS");
    YAML::Node params = YAML::LoadFile(workspaceDir+"/src/cito/config/params.yaml");
    std::string modelPathStr = workspaceDir + "/src/cito/model/" + params["model"].as<std::string>();
    const char *modelPath = modelPathStr.c_str();
    std::cout << "\n\nModel path: " << modelPath << "\n\n\n";
    // Load the model
    m = mj_loadXML(modelPath, NULL, NULL, 0);
    if( !m ) { mju_error("Cannot load the model"); }
    // Create and initialize data at the key pose
    d = mj_makeData(m);
    mju_copy(d->qpos, m->key_qpos, m->nq);
    // Evaluate the forward dynamics
    mj_forward(m, d);
    // Parse parameters from the model and the params file
    Params cp(m);

    // Create collision geometries and objects from the model
    std::unordered_map<int, fcl::CollisionObjectd*> fclObjects;
    for( int i=0; i<cp.nPair; i++ )
    {
        // Get the geometric type of the site
        for( int j=0; j<2; j++) {
            // Get the pose of the site
            fcl::Transform3d tf = getSiteTransform(d, cp.sites[i][j]);
            // Create the corresponding collision geometry and object
            std::shared_ptr<fcl::CollisionGeometryd> geom = createCollGeom(m, cp.sites[i][j], tf);
            fclObjects[cp.sites[i][j]] = new fcl::CollisionObjectd(geom, tf);
        }
    }
    // Compute distance
    for( int i=0; i<cp.nPair; i++ )
    {
        distRes.clear();
        fcl::distance(fclObjects[cp.sites[i][0]], fclObjects[cp.sites[i][1]], distReq, distRes);
        printf("Distance between sites %d and %d: %f.\n", cp.sites[i][0], cp.sites[i][1], distRes.min_distance);
        std::cout << "Nearest point on O1: " << distRes.nearest_points[0].transpose() << "\n";
        std::cout << "Nearest point on O2: " << distRes.nearest_points[1].transpose() << "\n";
    }

    // FCL test with random shapes
    std::cout << "\n\nCreate random shapes and calculate distance:\n";
    // Create shapes and compute their volumes
    std::shared_ptr<fcl::ShapeBased> shape1 = std::make_shared<fcl::Boxd>(1., 2., 3.);
    std::shared_ptr<fcl::ShapeBased> shape2 = std::make_shared<fcl::Sphered>(1.);
    // Object poses
    fcl::Transform3d tfS1, tfS2, tfS3;
    tfS1.translation() = fcl::Vector3d(5., 5., 5.);
    tfS2.translation() = fcl::Vector3d(0., 0., 0.);
    tfS3.translation() = fcl::Vector3d(10., 0., 0.);
    tfS1.linear().setIdentity();
    tfS2.linear().setIdentity();
    tfS3.linear().setIdentity();
    std::vector<fcl::Transform3d> transforms;
    transforms.push_back(tfS1);
    transforms.push_back(tfS2);
    transforms.push_back(tfS3);
    // Create collision geometries
    std::vector<std::shared_ptr<fcl::CollisionGeometryd>> coll_geoms;
    coll_geoms.push_back(shape1);
    coll_geoms.push_back(shape2);
    coll_geoms.push_back(std::make_shared<fcl::Boxd>(1., 1., 1.));
    // Create collision objects
    std::vector<fcl::CollisionObjectd*> coll_objs;
    printf("  Collision geometries:\n");
    for(int i=0; i<coll_geoms.size(); i++)
    {
        coll_objs.push_back(new fcl::CollisionObjectd(coll_geoms[i], transforms[i]));
        // Print shape properties
        std::cout << "\tShape " << i << " - node type: " << coll_geoms[i]->getNodeType() <<
                     ", volume: " << coll_geoms[i]->computeVolume() <<
                     ", pos: " << transforms[i].translation().transpose() << "\n";
    }
    // Modify the pose of a collision object
    fcl::Transform3d tfNew;
    tfNew.translation() = fcl::Vector3d(1., 2., 3.); tfNew.linear().setIdentity();
    coll_objs[1]->setTransform(tfNew);
    // Print collision object data
    printf("  Collision objects:\n");
    for(auto it : coll_objs)
    {
        std::cout << "\tPos: " << it->getTranslation().transpose() <<
                     ", object type: " << it->getNodeType() << "\n";
    }
    // Compute distance
    for(int i=0; i<coll_objs.size(); i++)
    {
        for(int j=0; j<coll_objs.size(); j++)
        {
            if(i!=j)
            {
                distRes.clear();
                fcl::distance(coll_objs[i], coll_objs[j], distReq, distRes);
                // Print results
                std::cout << "  Objects " << i << " & " << j <<
                            "\n\tDistance: " << distRes.min_distance <<
                            "\n\tNearest point on O" << i << ": " << distRes.nearest_points[0].transpose() <<
                            "\n\tNearest point on O" << j << ": " << distRes.nearest_points[1].transpose() << "\n";
            }
        }
    }

    // Clean up memory
    for(auto it : fclObjects)
        delete it.second;
    for(int i=0; i<coll_objs.size(); i++)
        delete coll_objs[i];
    mj_deleteData(d);
    mj_deleteModel(m);
    mj_deactivate();
}