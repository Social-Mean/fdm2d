# cmake_minimum_required(VERSION 3.10)
# project(fdm2d)
# set(FDM2D_SRC ${CMAKE_CURRENT_SOURCE_DIR}/fdm2d.cpp)
add_library(fdm2d fdm2d.cpp)
find_package(Eigen3 REQUIRED)
target_link_libraries(fdm2d PUBLIC Eigen3::Eigen stdc++fs)
# 在 thirdpart/eigen-3.4.0 子目录的 CMakeLists.txt 中
# find_package(Eigen3 REQUIRED)