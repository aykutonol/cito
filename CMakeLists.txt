cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)

set(MUJOCO_DIR "$ENV{MJ_HOME}")
set(SNOPT_DIR  "$ENV{SN_HOME}")

MESSAGE(SNOPT_DIR ${SNOPT_DIR})
###################################
## catkin specific configuration ##rosrun 
###################################
catkin_package(
 INCLUDE_DIRS include ${EIGEN3_INCLUDE_DIR}
)

###########
## Build ##
###########
## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${MUJOCO_DIR}/include
 ${SNOPT_DIR}/include
)

## Declare a C++ library
add_library(${PROJECT_NAME}
 	src/cito_control.cpp
	src/cito_numdiff.cpp
	src/cito_scvx.cpp
	src/cito_sqopt.cpp
	src/mj_savelog.cpp
)

################
## Node names ##
################
set(NODE_NAME1 main)
set(NODE_NAME2 playlog)

#################
## Executables ##
#################
# main: runs successive convexification w/ SQOPT and numerical differentiation in MuJoCo
add_executable(${NODE_NAME1} src/main.cpp)
### mujoco w/o rendering
target_link_libraries(${NODE_NAME1} ${MUJOCO_DIR}/bin/libmujoco200nogl.so
				    pthread
				    ${PROJECT_NAME}
				    ${SNOPT_DIR}/lib/libsnopt7_cpp.so
				    yaml-cpp
)
# playlog: shows motions generated
add_executable(${NODE_NAME2} src/mj_playlog.cpp)
### mujoco w/ rendering
target_link_libraries(${NODE_NAME2} ${MUJOCO_DIR}/bin/libmujoco200.so
				    ${MUJOCO_DIR}/bin/libglew.so
				    ${MUJOCO_DIR}/bin/libglfw.so.3
				    GL
				    yaml-cpp
)

