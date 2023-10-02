#include <fdm2d/fdm2d.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace {
namespace py = pybind11;

PYBIND11_MODULE(fdm, m) {
  py::class_<FDM2D>(m, "FDM2D")
      .def(py::init<int, int>())
      .def("solve", &FDM2D::solve)
      .def("getPhi", &FDM2D::getPhi)
      .def("add_D_BC", &FDM2D::add_D_BC)
      .def("add_D_BCs", &FDM2D::add_D_BCs);
}
}  // namespace