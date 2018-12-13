cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)

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
 include/mujoco
 ${EIGEN3_INCLUDE_DIR}
)

## Declare a C++ library
add_library(${PROJECT_NAME}
 	src/cito_control.cpp
	src/cito_numdiff.cpp
)
#add_library(cito_cont src/cito_cont.cpp)

################
## Node names ##
################
set(NODE_NAME1 flymanoid_test)
set(NODE_NAME2 flymanoid)

#################
## Executables ##
#################
add_executable(${NODE_NAME1} src/flymanoid_test.cpp)
# mujoco w/ rendering
target_link_libraries(${NODE_NAME1} ${PROJECT_SOURCE_DIR}/bin/libmujoco200.so ${PROJECT_SOURCE_DIR}/bin/libglew.so ${PROJECT_SOURCE_DIR}/bin/libglfw.so.3 GL)


add_executable(${NODE_NAME2} src/flymanoid.cpp)
# mujoco w/o rendering
target_link_libraries(${NODE_NAME2} ${PROJECT_SOURCE_DIR}/bin/libmujoco200nogl.so
				    pthread
				    ${PROJECT_NAME}
)

