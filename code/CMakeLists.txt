cmake_minimum_required(VERSION 2.8.12)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

FIND_PACKAGE(PythonInterp 3)

project(crisp CXX)
add_subdirectory(pybind11)
pybind11_add_module(crisp python.cpp lbp.cpp gibbs.cpp distribution.cpp crisp.cpp node.cpp factor.cpp seir.cpp)

add_executable (testit main.cpp node.cpp factor.cpp seir.cpp distribution.cpp)
