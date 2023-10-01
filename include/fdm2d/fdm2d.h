#pragma once
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/gmres.hpp>
#include <array>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// #define AMGCL_NO_BOOST
AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND()
struct TupleHash {
  template <typename T, typename U>
  std::size_t operator()(const std::tuple<T, U>& tuple) const {
    auto hash1 = std::hash<T>{}(std::get<0>(tuple));
    auto hash2 = std::hash<U>{}(std::get<1>(tuple));
    return hash1 ^ hash2;
  }
};

class FDM2D {
 public:
  FDM2D() = delete;
  FDM2D(int Nx, int Ny);
  void solve();
  Eigen::MatrixXd getPhi() const;

 private:
  int Nx;
  int Ny;
  int npts;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  Eigen::VectorXd x;
  Eigen::MatrixXd phi;
  std::unordered_map<int, std::tuple<int, int> > tag2idx;
  std::unordered_map<std::tuple<int, int>, int, TupleHash> idx2tag;

  void initMaps();
  void setMatrix();
  void setPoisson();
  void setBCs();
  void set_D_BCs();
  void set_N_BCs();
  void set_R_BCs();
  void solveMatrix();
  void postProcess();
  void eraseRow(Eigen::SparseMatrix<double>&, int);
  void set_D_BC(int i, int j, double value);
  void set_D_BC(int tag, double value);
  void set_D_BC(std::tuple<int, int> idx, double value);
  void saveMatrix(const Eigen::MatrixXd matrix, std::string filename) const;
};