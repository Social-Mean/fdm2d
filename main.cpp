#include <fdm2d/fdm2d.h>

#include <fstream>
// #include <format>
#include <iostream>
#include <memory>

int main() {
  int Nx = 21;
  int Ny = 21;
  FDM2D fdm{Nx, Ny};
  fdm.solve();
  return 0;
}