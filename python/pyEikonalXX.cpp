#include "eikonalxx/version.hpp"
#include "include/pyGeometry2d.hpp"
#include "include/pyGeometry3d.hpp"
#include "eikonalxx/enums.hpp"
#include <pybind11/pybind11.h>

PYBIND11_MODULE(pyEikonalXX, m)
{
    m.attr("__version__") = EIKONALXX_VERSION;
    m.attr("__doc__") = "A toolkit for solving the eikonal equation in seismic applications.";

    PEikonalXX::initializeGeometry2D(m);
    PEikonalXX::initializeGeometry3D(m);

    pybind11::enum_<EikonalXX::Ordering2D> (m, "Ordering2D")
        .value("natural", EikonalXX::Ordering2D::NATURAL,
               "The natural grid ordering in 2D.")
        .value("zx", EikonalXX::Ordering2D::ZX,
               "The natural grid ordering in 2D where x is the fastest changing dimension and z the slowest changing dimension.")
        .value("xz", EikonalXX::Ordering2D::XZ,
               "The transpose grid ordering in 2D where x is the slowest changing dimension and z the fastest changing dimension.");

    pybind11::enum_<EikonalXX::Ordering3D> (m, "Ordering3D")
        .value("natural", EikonalXX::Ordering3D::NATURAL,
               "The natural grid ordering in 3D.")
        .value("zyx", EikonalXX::Ordering3D::ZYX,
               "The natural grid ordering in 3D where x is the fastest changing dimension, y the intermediate dimension, and z the slowest changing dimension.")
        .value("xyz", EikonalXX::Ordering3D::XYZ,
               "The transpose grid ordering in 3D where x is the slowest changing dimension, y the intermediate dimension, and z the fastest changing dimension.");


    pybind11::enum_<EikonalXX::SolverAlgorithm> (m, "SolverAlgorithm")
        .value("level_set_method", EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD,
               "The solver solves the eikonal equation using fast-sweeping where nodes in a level-set can be updated in parallel.  The parallelism can result in obtaining a solution faster but with considerably higher memory overhead.")
        .value("fast_sweeping_method", EikonalXX::SolverAlgorithm::FAST_SWEEPING_METHOD,
               "The solver solves the eikonal equation using fast-sweeping.  All updates during the sweep are performed serially.");

}
