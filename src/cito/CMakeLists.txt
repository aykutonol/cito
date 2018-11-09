cmake_minimum_required(VERSION 2.8.3)
project(cito)

## Compile as C++11, supported in ROS Kinetic and newer
 add_compile_options(-std=c++11)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED)

################################################
## Declare ROS dynamic reconfigure parameters ##
################################################

###################################
## catkin specific configuration ##
###################################
catkin_package(
 
)

###########
## Build ##
###########
## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
 include
# ${catkin_INCLUDE_DIRS}
)

## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/cito.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})


################
## Node names ##
################
set(NODE_NAME1 flymanoid_test)

#################
## Executables ##
#################
add_executable(${NODE_NAME1} src/flymanoid_test.cpp)
target_link_libraries(${NODE_NAME1} ${PROJECT_SOURCE_DIR}/bin/libmujoco200.so ${PROJECT_SOURCE_DIR}/bin/libglew.so ${PROJECT_SOURCE_DIR}/bin/libglfw.so.3 GL)


