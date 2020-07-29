cmake_minimum_required(VERSION 2.8)
project(oed)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -lnlopt -lm")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wno-unused -Wno-unused-parameter -Wno-unknown-pragmas -Wno-ignored-attributes")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef -Werror -Wno-unused -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors")

INCLUDE_DIRECTORIES ( "src" "src/merr" "src/pce" "extern/eigen" )

set(SHARED_SOURCES
    src/RandomGenerator.cpp
    src/MultivariateNormal.cpp
    src/MultivariateT.cpp
    src/Utilities.cpp)

set(EXPDESIGN_SOURCES
    src/ExpectedInformationEstimator.cpp
    src/LinearGaussianExperiment.cpp)

add_executable(Driver src/Driver.cpp ${SHARED_SOURCES} ${EXPDESIGN_SOURCES})
# add_executable(DriverMC src/DriverMC.cpp ${SHARED_SOURCES} ${EXPDESIGN_SOURCES})

set(MERR_SOURCES
    src/merr/Model.cpp
    src/merr/ModelInference.cpp)

# add_executable(ModelInference ${MERR_SOURCES} ${SHARED_SOURCES})
