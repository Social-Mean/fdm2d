#include <fdm2d/fdm2d.h>

FDM2D::FDM2D(int Nx_, int Ny_) : Nx{Nx_}, Ny{Ny_}, npts{Nx_ * Ny_} {
  A.resize(npts, npts);
  b.resize(npts);
  x.resize(npts);
  phi.resize(Nx, Ny);
  initMaps();
}

void FDM2D::initMaps() {
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      int tag = i * Ny + j;
      std::tuple idx{i, j};
      tag2idx[tag] = idx;
      idx2tag[idx] = tag;
    }
  }
}

void FDM2D::setMatrix() {
  setPoisson();
  setBCs();
}

void FDM2D::setPoisson() {
  using T = Eigen::Triplet<double>;
  int listSize = 5 * (Nx - 2) * (Ny - 2);
  std::vector<T> tripletList;
  tripletList.reserve(listSize);
  for (int i = 1; i < Nx - 1; ++i) {
    for (int j = 1; j < Ny - 1; ++j) {
      auto idx{std::make_tuple(i, j)};
      auto tag = idx2tag[idx];
      tripletList.push_back(T(tag, tag, -4));
      tripletList.push_back(T(tag, idx2tag[std::make_tuple(i - 1, j)], 1));
      tripletList.push_back(T(tag, idx2tag[std::make_tuple(i + 1, j)], 1));
      tripletList.push_back(T(tag, idx2tag[std::make_tuple(i, j - 1)], 1));
      tripletList.push_back(T(tag, idx2tag[std::make_tuple(i, j + 1)], 1));
      b(tag) = 0;
    }
  }
  A.setFromTriplets(tripletList.begin(), tripletList.end());
}

void FDM2D::setBCs() {
  set_D_BCs();
  set_N_BCs();
  set_R_BCs();
}

void FDM2D::set_D_BCs() {
  // for (int i = 0; i < Nx; ++i) {
  //   for (int j : {0, Ny - 1}) {
  //     set_D_BC(i, j, 0);
  //   }
  // }

  // for (int j = 0; j < Ny; ++j) {
  //   for (int i : {0, Nx - 1}) {
  //     set_D_BC(i, j, 0);
  //   }
  // }

  // set_D_BC(Nx / 2, Ny / 2, 1);
  for (auto& [i, j, value] : D_BCs) {
    set_D_BC(i, j, value);
  }

  A.makeCompressed();
}

void FDM2D::set_N_BCs() {}

void FDM2D::set_R_BCs() {}

void FDM2D::solve() {
  setMatrix();
  solveMatrix();
  postProcess();
}

void FDM2D::solveMatrix() {
  amgcl::profiler<> prof;
  prof.tic("solve");
  std::cout << "开始求解" << std::endl;
  // {  // 使用 amgcl 的求解器
  //   using Solver = amgcl::make_solver<
  //       amgcl::amg<amgcl::backend::builtin<double>,
  //                  amgcl::coarsening::smoothed_aggregation,
  //                  amgcl::relaxation::spai0>,
  //       amgcl::solver::bicgstab<amgcl::backend::builtin<double> > >;
  //   Solver solver(A);
  //   auto [iter, error] = solver(b, x);
  //   std::cout << iter << std::endl << error << std::endl;
  // }
  {  // 使用 Eigen 自带的求解器
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    solver.setTolerance(1e-5);
    solver.compute(A);
    x = solver.solve(b);
    if (solver.info() == Eigen::Success) {
      std::cout << "Complete Solving!" << std::endl;
    }
  }

  phi = x.reshaped<Eigen::RowMajor>(phi.rows(), phi.cols());
  prof.toc("solve");

  std::cout << prof << std::endl;
}

void FDM2D::postProcess() {}

Eigen::MatrixXd FDM2D::getPhi() const { return phi; }

void FDM2D::eraseRow(Eigen::SparseMatrix<double>& matrix, int _row) {
  matrix.prune([_row](int row, int col, double value) { return row != _row; });
}

void FDM2D::set_D_BC(int i, int j, double value) {
  auto tag = idx2tag[std::make_tuple(i, j)];
  eraseRow(A, tag);
  A.insert(tag, tag) = 1;
  b(tag) = value;
}

void FDM2D::saveMatrix(const Eigen::MatrixXd matrix,
                       std::string filename) const {
  std::cout << std::format("开始保存\"{}\"\t", filename);
  std::cout.flush();
  if (std::ofstream file{filename}; file.is_open()) {
    file << matrix;
  }
  std::cout << "保存完成" << std::endl;
}

void FDM2D::add_D_BC(int i, int j, double value) {
  D_BCs.push_back(std::make_tuple(i, j, value));
}

void FDM2D::add_D_BCs(const std::vector<int>& is, const std::vector<int>& js,
                      const std::vector<double>& values) {
  // if (!(is.size() == js.size() && is.size() == values.size())) {
  // // BUG: pybind11 的异常处理
  //   throw "数组长度不匹配";
  // }
  // TODO: 创建一个向量, 存储此次增加的边界条件, 再连接D_BCs和这个新增加的向量
  for (auto [i, j, value] : std::views::zip(is, js, values)) {
    add_D_BC(i, j, value);
  }
  // for (auto i = std::begin)
}