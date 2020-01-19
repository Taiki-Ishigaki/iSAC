add_library(sac
  SHARED
  sac.cpp
  )

target_include_directories(sac
  PUBLIC ${PROJECT_SOURCE_DIR}/include
  )