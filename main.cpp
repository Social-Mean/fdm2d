#include <fdm2d/fdm2d.h>

#include <fstream>
// #include <format>
#include <iostream>
#include <memory>

int main() {
  int Nx = 11;
  int Ny = 11;
  std::cin >> Nx >> Ny;
  auto fdm = std::make_unique<FDM2D>(Nx, Ny);
  // FDM2D fdm{Nx, Ny};
  fdm->solve();
  return 0;
}