add_library(sac
  SHARED
  sac.cpp
  )

target_include_directories(sac
  PUBLIC ${PROJECT_SOURCE_DIR}/include
  )

find_package(Eigen3 REQUIRED)
MESSAGE( [Main] " EIGEN3_INCLUDE_DIRS = ${EIGEN3_INCLUDE_DIRS}") 
# EIGEN3_INCLUDE_DIRS = /usr/include/eigen3

include_directories(
    ${EIGEN3_INCLUDE_DIRS}
)