# Contact-Implicit Trajectory Optimization (CITO)

This package is developed for planning non-prehensile manipulation and locomotion
motions without a predefined contact schedule. The optimization-based planning algorithm 
is based on a variable smooth contact model and successive convexification. 
Please see [1] for a detailed description of the algorithm.

In this framework, MuJoCo is used to evaluate the nonlinear dynamics. The partial derivatives
of the dynamics about the previous trajectory are obtained by numerical differentiation.
The resulting convex subproblems are solved by SQOPT since the dynamic constraints and
the cost are sparse.

The libraries are implemented in C++ and catkin is used to compile them.

**Please note that both the method and the code are currently under development.**

## Dependencies
- [ROS (catkin)](http://wiki.ros.org/catkin)
- [Eigen 3](https://eigen.tuxfamily.org/dox/GettingStarted.html)
- [MuJoCo 2.00](http://www.mujoco.org/)
- [SNOPT 7](https://ccom.ucsd.edu/~optimizers/solvers/snopt/)
- [YAML](https://github.com/jbeder/yaml-cpp)

## Installation
1. Create a workspace and download the code:  
    ```
    mkdir -p ~/cito_ws/src
    cd ~/cito_ws/src/
    git clone https://github.com/aykutonol/cito.git
    ```  
2. Set the environment variables:  
    ```
    export CITO_WS=~/cito_ws  
    export MJ_KEY=*path to the licence file for MuJoCo*
    export MJ_HOME=*path to the home directory of MuJoCo*  
    export SN_HOME=*path to the home directory of SNOPT*  
    ```
3. Build the package:
    ```
    cd ~/cito_ws/
    catkin build
    source devel/setup.bash
    ```

## Usage
The model file and the application-specific parameters related to the model (i.e., joint types,
contact pairs, etc.) and simulation (i.e., time horizon length and dynamic and control time step
sizes) are  defined in include/cito_params.h.

A motion can be planned by running the main:  
`rosrun cito main`  
which will record the optimal trajectory into cito_ws/logs/mjLog_model_name, a binary file.

The planned motion can be played back by:  
`rosrun cito playlog`

While rendering, the trajectories (i.e., time, joint positions, velocities, accelerations, and 
actuation forces) are recorded into logs/traj_model_name (i.e., a readable file) so that they can 
be sent to the robot.

In order to change the model, include/cito_params.h and config/task.yaml need to be modified.

The task can be modified through config/task.yaml.

The successive convexification parameters and the variable smooth contact model parameters can be 
changed through config/scvx.yaml and config/vscm.yaml.


## Examples
**Sawyer tabletop pushing**  
In this application, the task is to push a box on a table with a 7 DOF robot arm, Sawyer.
The model for this example is model/sawyer_push.xml. To see an example motion, you need to:
- copy include/cito_params_sawyer.h and paste into include/cito_params,
- copy config/task_sawyer.yaml and paste into config/task.yaml.

**Simple humanoid locomotion at zero gravity (_Flymanoid_)**  
The goal of this application is to plan a locomotion behavior for a planar human-like robot.
The model for this example is model/flymanoid.xml. To see an example motion, you need to:  
- copy include/cito_params_flymanoid.h and paste into include/cito_params,
- copy config/task_flymanoid.yaml and paste into config/task.yaml.


## Citing
If you use this package, please cite the following paper:

[1] [\"{O}nol, A. \"{O}., Long, P., & Padir, T. (2019). Contact-Implicit TrajectoryOptimization
Based on a Variable Smooth Contact Model and Successive Convexification.
In _2019 IEEE International Conference on Robotics and Automation (ICRA)_. IEEE.](https://arxiv.org/abs/1810.10462
) [Accepted]
```
@inproceedings{onol2019contact,
  title={Contact-Implicit Trajectory Optimization Based on a Variable Smooth Contact Model and Successive Convexification},
  author={Onol, Aykut Ozgun and Long, Philip and Padir, Taskin},
  booktitle={2019 IEEE International Conference on Robotics and Automation (ICRA)},
  year={2019}
}
```

## TODOS
- Integrate FCL for distance calculation
- Analytic derivatives
