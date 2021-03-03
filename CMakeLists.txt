cmake_minimum_required(VERSION 3.13)
project(EMPSIPON C CXX)

set(lib_source_files
    src/Auxiliary.cpp
    src/b3D_integrator.cpp
    src/distributions.cpp
    src/domain.cpp
    src/dump.cpp
    src/grid.cpp
    src/hopping.cpp
    src/netmin.cpp
    src/network.cpp
    src/non_bonded_scheme_routines.cpp
    src/rng.cpp)

add_library(EmpsiponCore ${lib_source_files})
target_include_directories(EmpsiponCore PUBLIC "src/include")

add_executable(netmin.x src/main.cpp)
target_link_libraries(netmin.x PRIVATE EmpsiponCore)

file(COPY ./configurations
    DESTINATION ${CMAKE_BINARY_DIR})
