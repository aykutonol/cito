cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp 0.6 REQUIRED)

set(MUJOCO_DIR    "$ENV{MJ_HOME}")
set(SNOPT_DIR     "$ENV{SN_HOME}")

# Avoid CPM0042 warning on OSX
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CMAKE_MACOSX_RPATH FALSE)
endif()

catkin_package()

include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${MUJOCO_DIR}/include
 ${SNOPT_DIR}/include
)

# Create the CITO library
add_library(${PROJECT_NAME}
	src/cito_params.cpp
	src/cito_control.cpp
	src/cito_numdiff.cpp
	src/cito_scvx.cpp
	src/cito_sqopt.cpp
	src/mj_savelog.cpp
)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                          pthread
                          yaml-cpp
	)
endif()

# Set the executable names
set(NODE_NAME1 main)
set(NODE_NAME2 playlog)

# Create executables
# main: solves the CITO problem using successive convexification
#       w/ SQOPT and numerical differentiation in MuJoCo
add_executable(${NODE_NAME1} src/main.cpp)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${NODE_NAME1}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${NODE_NAME1}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                          pthread
                          yaml-cpp
	)
endif()

# playlog: renders the generated motion
add_executable(${NODE_NAME2} src/mj_playlog.cpp)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	target_link_libraries(${NODE_NAME2}
                          ${MUJOCO_DIR}/bin/libmujoco200.so
                          ${MUJOCO_DIR}/bin/libglew.so
                          ${MUJOCO_DIR}/bin/libglfw.so.3
                          GL
                          pthread
                          yaml-cpp
	)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	target_link_libraries(${NODE_NAME2}
                          ${MUJOCO_DIR}/bin/libmujoco200.dylib
                          ${MUJOCO_DIR}/bin/libglfw.3.dylib
                          pthread
                          yaml-cpp
	)
endif()