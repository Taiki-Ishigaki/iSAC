cmake_minimum_required(VERSION 2.8)
project(tank_controller CXX)

# find_package(PkgConfig REQUIRED)
# pkg_check_modules(CHOREONOID choreonoid-body REQUIRED)

find_package( Eigen3 REQUIRED)
MESSAGE( [Main] " EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}") 
# EIGEN3_INCLUDE_DIRS = /usr/include/eigen3

find_package( choreonoid REQUIRED)
MESSAGE( [Main] " CHOREONOID_INCLUDE_DIRS = ${CHOREONOID_INCLUDE_DIRS}") 

include_directories(
    ${EIGEN3_INCLUDE_DIRS}
    ${CHOREONOID_INCLUDE_DIRS}
)

link_directories(
    ${CHOREONOID_LIBRARY_DIRS}
)

add_library(tank_controller
    SHARED
    sac.cpp
    isac.cpp
    tank_controller.cpp
)

target_link_libraries(tank_controller
    ${CHOREONOID_BODY_LIBRARIES}
)

set_target_properties(tank_controller PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CHOREONOID_PLUGIN_SUBDIR}/simplecontroller PREFIX "")

install(TARGETS tank_controller LIBRARY DESTINATION ${CHOREONOID_PLUGIN_SUBDIR}/simplecontroller)