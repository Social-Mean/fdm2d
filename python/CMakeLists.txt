cmake_minimum_required(VERSION 3.5...3.26)

pybind11_add_module(fdm fdm2d_wrapper.cpp)
target_link_libraries(fdm PUBLIC fdm2d Eigen3::Eigen amgcl::amgcl eigen_target stdc++fs)
