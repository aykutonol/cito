cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)

set(MUJOCO_DIR    "$ENV{MJ_HOME}")
set(SNOPT_DIR     "$ENV{SN_HOME}")
set(PINOCCHIO_DIR "$ENV{PN_HOME}")

catkin_package()

include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${MUJOCO_DIR}/include
 ${SNOPT_DIR}/include
 ${PINOCCHIO_DIR}/build/include
)

add_library(${PROJECT_NAME}
 	src/cito_params.cpp
	src/cito_control.cpp
	src/cito_numdiff.cpp
	src/cito_scvx.cpp
	src/cito_sqopt.cpp
	src/mj_savelog.cpp
)

set(NODE_NAME1 main)
set(NODE_NAME2 playlog)
set(NODE_NAME3 deriv)

# main: runs successive convexification w/ SQOPT and numerical differentiation in MuJoCo
add_executable(${NODE_NAME1} src/main.cpp)
### MuJoCo w/o rendering
target_link_libraries(${NODE_NAME1} ${MUJOCO_DIR}/bin/libmujoco200nogl.so
				    pthread
				    ${PROJECT_NAME}
				    ${SNOPT_DIR}/lib/libsnopt7_cpp.so
				    yaml-cpp
)
# playlog: shows motions generated
add_executable(${NODE_NAME2} src/mj_playlog.cpp)
### MuJoCo w/ rendering
target_link_libraries(${NODE_NAME2} ${MUJOCO_DIR}/bin/libmujoco200.so
				    ${MUJOCO_DIR}/bin/libglew.so
				    ${MUJOCO_DIR}/bin/libglfw.so.3
				    GL
				    yaml-cpp
)
# deriv: compares Pinocchio and MuJoCo derivatives
add_executable(${NODE_NAME3} src/derivative_comparison.cpp)
target_link_libraries(${NODE_NAME3} ${MUJOCO_DIR}/bin/libmujoco200nogl.so
					pthread
					${PROJECT_NAME}
					yaml-cpp
					${SNOPT_DIR}/lib/libsnopt7_cpp.so
					${Boost_SYSTEM_LIBRARY}
					urdfdom_model
)
target_compile_definitions(${NODE_NAME3} PUBLIC PINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
										 	    PINOCCHIO_WITH_URDFDOM
)