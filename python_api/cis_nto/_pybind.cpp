// pybind11 bindings for the C++ EasyWFO core.
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/complex.h>
#include "easywfo_core.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "C++ core for cis_nto: Shakiba & Akimov reprojection matrix";
    m.def("orthogonalized_T", &EasyWFO::orthogonalized_T,
          py::arg("T"), py::arg("tol") = 1e-10, py::arg("drop") = false);
    m.def("truncate_wavefunction", &EasyWFO::truncate_wavefunction,
          py::arg("coeffs"), py::arg("wfthres") = 0.99);
    m.def("propagate_reprojection", &EasyWFO::propagate_reprojection,
          py::arg("C_old"), py::arg("H"), py::arg("S_MCH"), py::arg("dt"));
}
