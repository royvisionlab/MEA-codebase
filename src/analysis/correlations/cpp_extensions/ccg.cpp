#include <vector>
#include <cstdint>
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "ccg.hpp"

namespace py = pybind11;

PYBIND11_MODULE(correlation_cpp_extensions, m) {

    m.doc() = "pybind11 module for ccg.cpp"; // optional module docstring

    m.def("correlation_coefficient",
            &correlation_coefficient,
            "Compute the correlation coefficient, ignoring NaNs.");

    m.def("crosscorrelogram",
            &crosscorrelogram,
            pybind11::return_value_policy::take_ownership,
            "Compute spike time correlations between two spike trains.");

    m.def("autocorrelogram",
            &autocorrelogram,
            pybind11::return_value_policy::take_ownership,
            "Compute spike time autocorrelogram for a spike train.");
}
