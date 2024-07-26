cmake_minimum_required(VERSION 3.11)
add_definitions(-std=c++11)

project(noorbs)
set(CMAKE_MODULE_PATH ${Python_INCLUDE_DIRS}/cmake)

# options
option(BUILD_PYTHON_BINDINGS "create python bindings" ON)

## eigen
find_package(Eigen3)
include_directories(${EIGEN3_INCLUDE_DIR})

if(BUILD_PYTHON_BINDINGS)
    find_package(Python COMPONENTS Interpreter Development)
    find_package(pybind11 CONFIG)
endif()

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
endif()

add_subdirectory("src/")

if (BUILD_TESTS)
    add_subdirectory("tests/")
endif(BUILD_TESTS)
