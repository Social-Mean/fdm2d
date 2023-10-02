cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(BOOST_INCLUDE_DIR /usr/include)
set(BOOST_ROOT /usr/lib64)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

project(fdm LANGUAGES CXX)

add_library(eigen_target INTERFACE)

include_directories(./include)
include_directories(./thirdpart/amgcl)

add_subdirectory(./src)
add_subdirectory(./thirdpart/eigen)
add_subdirectory(./thirdpart/amgcl)
add_subdirectory(./thirdpart/pybind11)
add_subdirectory(./python)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)