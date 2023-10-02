#include <fdm2d/fdm2d.h>
#include <pybind11/pybind11.h>

namespace {
namespace py = pybind11;

PYBIND11_MODULE(fdm, m) {
  py::class_<FDM2D>(m, "FDM2D")
      .def(py::init<int, int>())
      .def("solve", &FDM2D::solve)
      .def("getPhi", &FDM2D::getPhi);
}
}  // namespace