add_executable(test.out main.cpp)

find_package(Eigen3 REQUIRED)
MESSAGE( [Main] " EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}") 
# EIGEN3_INCLUDE_DIRS = /usr/include/eigen3

include_directories(
    ${EIGEN3_INCLUDE_DIRS}
)

target_link_libraries(test.out sac)