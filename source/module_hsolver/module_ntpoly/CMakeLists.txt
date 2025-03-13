add_library(ntpoly OBJECT simple_ntpoly.cpp utils.cpp)

if(ENABLE_COVERAGE)
  add_coverage(ntpoly)
endif()
