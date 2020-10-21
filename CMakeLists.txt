cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(fcl REQUIRED)

# Get environment variables
set(MUJOCO_DIR    "$ENV{MJ_HOME}")
set(SNOPT_DIR     "$ENV{SN_HOME}")

# Avoid CPM0042 warning on OSX
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CMAKE_MACOSX_RPATH FALSE)
endif()

# Set build type and compilation flags
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic")

catkin_package()

include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${MUJOCO_DIR}/include
 ${SNOPT_DIR}/include
)

# Create the CITO library
add_library(${PROJECT_NAME}
	src/cito/params.cpp
	src/cito/control.cpp
	src/cito/numdiff.cpp
	src/cito/scvx.cpp
	src/cito/sqopt.cpp
	src/cito/savelog.cpp
	src/cito/penalty_loop.cpp
)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
                          fcl
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                          pthread
                          yaml-cpp
                          fcl
	)
endif()

# Set the executable names
set(NODE_NAME1 main)
set(NODE_NAME2 playlog)
set(NODE_NAME3 fcl_funcs)

# Create executables
# main: solves the CITO problem using successive convexification
#       w/ SQOPT and numerical differentiation in MuJoCo
add_executable(${NODE_NAME1} src/main.cpp)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${NODE_NAME1}
                          ${PROJECT_NAME}
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${NODE_NAME1}
                          ${PROJECT_NAME}
	)
endif()

# playlog: renders the generated motion
add_executable(${NODE_NAME2} src/playlog.cpp)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${NODE_NAME2}
                          ${MUJOCO_DIR}/bin/libmujoco200.so
                          ${MUJOCO_DIR}/bin/libglew.so
                          ${MUJOCO_DIR}/bin/libglfw.so.3
                          GL
                          yaml-cpp
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${NODE_NAME2}
                          ${MUJOCO_DIR}/bin/libmujoco200.dylib
                          ${MUJOCO_DIR}/bin/libglfw.3.dylib
                          yaml-cpp
	)
endif()

# fcl_funcs: develop and test fcl functions
add_executable(${NODE_NAME3} test/fcl_funcs.cpp)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${NODE_NAME3}
                          ${PROJECT_NAME}
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${NODE_NAME3}
                          ${PROJECT_NAME}
	)
endif()

# Unit tests
catkin_add_gtest(${PROJECT_NAME}-test test/test_cito.cpp)
if(TARGET ${PROJECT_NAME}-test)
    target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
endif()