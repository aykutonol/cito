cmake_minimum_required(VERSION 2.8.3)
project(cito)

add_compile_options(-std=c++11)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(yaml-cpp 0.6 REQUIRED)
find_package(pinocchio REQUIRED)
find_package(boost REQUIRED COMPONENTS system)

# On OSX, avoid CPM0042 warning and link /usr/local/lib
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH FALSE)
    link_directories("/usr/local/lib")
endif()

set(MUJOCO_DIR    "$ENV{MJ_HOME}")
set(SNOPT_DIR     "$ENV{SN_HOME}")

catkin_package()

include_directories(
 include
 ${EIGEN3_INCLUDE_DIR}
 ${MUJOCO_DIR}/include
 ${SNOPT_DIR}/include
 ${pinocchio_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
     src/cito_params.cpp
    src/cito_control.cpp
    src/cito_numdiff.cpp
    src/cito_deriv.cpp
    src/cito_scvx.cpp
    src/cito_sqopt.cpp
    src/mj_savelog.cpp
)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
                          ${Boost_SYSTEM_LIBRARY}
                          urdfdom_model
    )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                          pthread
                          yaml-cpp
                          ${Boost_SYSTEM_LIBRARY}
                          urdfdom_model
    )
endif()
target_compile_definitions(${PROJECT_NAME} PUBLIC PINOCCHIO_URDFDOM_TYPEDEF_SHARED_PTR
                                                  PINOCCHIO_URDFDOM_USE_STD_SHARED_PTR
                                                  BOOST_MPL_LIMIT_LIST_SIZE=30
                                                  PINOCCHIO_WITH_URDFDOM
)

# Set executable names
set(NODE_NAME1 main)
set(NODE_NAME2 playlog)
set(NODE_NAME3 tosc)
set(NODE_NAME4 deriv)
set(NODE_NAME5 contact)

# Create executables
# main: runs SCVX w/ numerical differentiation in MuJoCo
add_executable(${NODE_NAME1} src/main.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
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

# playlog: renders generated motions
add_executable(${NODE_NAME2} src/mj_playlog.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
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

# tosc: runs SCVX w/ analytical derivatives from Pinocchio
add_executable(${NODE_NAME3} src/tosc.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(${NODE_NAME3}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
    )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${NODE_NAME3}
                            ${PROJECT_NAME}
                            ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                            ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                            pthread
                            yaml-cpp
    )
endif()

# deriv: compares Pinocchio and MuJoCo derivatives
add_executable(${NODE_NAME4} src/derivative_comparison.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(${NODE_NAME4}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
    )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${NODE_NAME4}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                          pthread
                          yaml-cpp
    )
endif()


# contact: compares Pinocchio and MuJoCo derivatives w/ contact
add_executable(${NODE_NAME5} src/derivative_contact.cpp)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    target_link_libraries(${NODE_NAME5}
                          ${PROJECT_NAME}
                          ${MUJOCO_DIR}/bin/libmujoco200nogl.so
                          ${SNOPT_DIR}/lib/libsnopt7_cpp.so
                          pthread
                          yaml-cpp
    )
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_link_libraries(${NODE_NAME5}
                            ${PROJECT_NAME}
                            ${MUJOCO_DIR}/bin/libmujoco200nogl.dylib
                            ${SNOPT_DIR}/lib/libsnopt7_cpp.dylib
                            pthread
                            yaml-cpp
    )
endif()