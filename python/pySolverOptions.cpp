#include <string>
#include <eikonalxx/solverOptions.hpp>
#include "include/pySolverOptions.hpp"

using namespace PEikonalXX;

/// C'tor
SolverOptions::SolverOptions() :
    pImpl(std::make_unique<EikonalXX::SolverOptions> ())
{
}

/// Copy c'tors
SolverOptions::SolverOptions(const SolverOptions &options)
{
    *this = options;
}

SolverOptions::SolverOptions(const EikonalXX::SolverOptions &options)
{
    *this = options;
}

/// Move c'tor
SolverOptions::SolverOptions(SolverOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Destructor
SolverOptions::~SolverOptions() = default;

/// Copy assignment
SolverOptions& SolverOptions::operator=(const SolverOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<EikonalXX::SolverOptions> (*options.pImpl);
    return *this;
}

/// Copy assignment
SolverOptions& SolverOptions::operator=(const EikonalXX::SolverOptions &options)
{
    pImpl = std::make_unique<EikonalXX::SolverOptions> (options);
    return *this;
}

/// Move assignment
SolverOptions& SolverOptions::operator=(SolverOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Reset class
void SolverOptions::clear() noexcept
{
    pImpl->clear();
}

/// Number of sweeps
void SolverOptions::setNumberOfSweeps(const uint16_t nSweeps) noexcept
{
    pImpl->setNumberOfSweeps(nSweeps);
}

int SolverOptions::getNumberOfSweeps() const noexcept
{
    return pImpl->getNumberOfSweeps();
}
 
/// Convergence tolerance
void SolverOptions::setConvergenceTolerance(const double tolerance) noexcept
{
    pImpl->setConvergenceTolerance(tolerance);
}

double SolverOptions::getConvergenceTolerance() const noexcept
{
    return pImpl->getConvergenceTolerance();
}

/// Spherical solver radius
void SolverOptions::setSphericalSolverRadius(const int epsilon) noexcept
{
    pImpl->setSphericalSolverRadius(epsilon);
}

int SolverOptions::getSphericalSolverRadius() const noexcept
{
    return pImpl->getSphericalSolverRadius();
}

/// Verbosity
void SolverOptions::setVerbosity(EikonalXX::Verbosity verbosity) noexcept
{
    pImpl->setVerbosity(verbosity);
}
    
EikonalXX::Verbosity SolverOptions::getVerbosity() const noexcept
{
    return pImpl->getVerbosity();
}

/// Sovler algorithm
void SolverOptions::setAlgorithm(
    const EikonalXX::SolverAlgorithm algorithm) noexcept
{
    pImpl->setAlgorithm(algorithm);
}

EikonalXX::SolverAlgorithm SolverOptions::getAlgorithm() const noexcept
{
    return pImpl->getAlgorithm();
}

void PEikonalXX::initializeSolverOptions(pybind11::module &m)
{
    pybind11::class_<PEikonalXX::SolverOptions> o(m, "SolverOptions");
    o.def(pybind11::init<> ());
    o.doc() = "This defines the solver options which are as follows:\n\nnumber_of_sweeps : The number of Gauss-Seidel iterations used in the fast-sweeping method.  Convergence can be terminated early if the change in the travel time solution is less than the tolerance.\n\ntolerance : If the maximum change in the travel time solution (measured in seconds) is less than this value then the Gauss-Seidel iterations will be terminated early.  This can be disabled by setting the value to zero or a negative number.\n\nspherical_solver_radius : When the update node's grid index in a principle direction exceeds this value then the solver will use a Cartesian finite difference stencil as opposed to a factored-eikonal equation finite-difference stencil designed to well-approximate the solution in the presence of high-wavefront curvature.\n\nalgorithm : The solver algorithm.\n\nverbosity : The solver's verbosity.";
    o.def_property("number_of_sweeps",
                   &SolverOptions::getNumberOfSweeps,
                   &SolverOptions::setNumberOfSweeps);
    o.def_property("tolerance",
                   &SolverOptions::getConvergenceTolerance,
                   &SolverOptions::setConvergenceTolerance);
    o.def_property("spherical_solver_radius",
                   &SolverOptions::getSphericalSolverRadius,
                   &SolverOptions::setSphericalSolverRadius);
    o.def_property("algorithm",
                   &SolverOptions::getAlgorithm,
                   &SolverOptions::setAlgorithm);
    o.def_property("verbosity",
                   &SolverOptions::getVerbosity,
                   &SolverOptions::setVerbosity);
    o.def("clear",
          &SolverOptions::clear,
          "Resets the class.");

    /// Pickling rules (makes this class copyable)
    o.def(pybind11::pickle(
        [](const SolverOptions &o) { //__getstate__
            auto nSweeps = static_cast<int> (o.getNumberOfSweeps());
            auto nEps = o.getSphericalSolverRadius();
            auto tol = o.getConvergenceTolerance();
            auto algorithm = static_cast<int> (o.getAlgorithm());
            auto verbosity = static_cast<int> (o.getVerbosity());
            return pybind11::make_tuple(tol, nSweeps, nEps,
                                        algorithm, verbosity);
        },
        [](pybind11::tuple t) { //__setstate__
            if (t.size() != 5)
            {
                throw std::runtime_error("Invalid state!");
            }
            auto tol = t[0].cast<double> (); 
            auto nSweeps = static_cast<uint16_t> (t[1].cast<int> ());
            auto nEps = t[2].cast<int> (); 
            auto algorithm
                = static_cast<EikonalXX::SolverAlgorithm> (t[3].cast<int> ());
            auto verbosity
                = static_cast<EikonalXX::Verbosity> (t[4].cast<int> ());
            SolverOptions p;
            p.setConvergenceTolerance(tol);
            p.setNumberOfSweeps(nSweeps);
            p.setSphericalSolverRadius(nEps);
            p.setAlgorithm(algorithm);
            p.setVerbosity(verbosity);
            return p;
        }
    ));

}
