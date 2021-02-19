#include "eikonalxx/solverOptions.hpp"

using namespace EikonalXX;

#define DEFFAULT_CONVERGENCE_TOLERANCE 1.e-5
#define DEFAULT_SPHERICAL_RADIUS 3
#define DEFAULT_NUMBER_OF_SWEEPS 5

class SolverOptions::SolverOptionsImpl
{
public:
    /// Convergence tolerance in seconds.
    double mConvergenceTolerance = DEFFAULT_CONVERGENCE_TOLERANCE;
    /// Spherical solver radius in grid points.
    int mSphericalRadius = DEFAULT_SPHERICAL_RADIUS;
    /// Solver algorithm
    SolverAlgorithm mAlgorithm = SolverAlgorithm::LEVEL_SET_METHOD;
    /// Verbosity
    Verbosity mVerbosity = Verbosity::ERROR;
    /// Number of sweeps (Gauss-Seidel iterations)
    uint16_t mSweeps = DEFAULT_NUMBER_OF_SWEEPS;
};

/// Constructor
SolverOptions::SolverOptions() :
    pImpl(std::make_unique<SolverOptionsImpl> ())
{
}

/// Copy c'tor
SolverOptions::SolverOptions(const SolverOptions &options)
{
    *this = options;
}

/// Solve c'tor
SolverOptions::SolverOptions(SolverOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Copy assignment
SolverOptions& SolverOptions::operator=(const SolverOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<SolverOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
SolverOptions& SolverOptions::operator=(SolverOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Destructor
SolverOptions::~SolverOptions() = default;

/// Sets the number of sweeps
void SolverOptions::setNumberOfSweeps(uint16_t nSweeps) noexcept
{
    pImpl->mSweeps = nSweeps;
}

/// Gets the number of sweeps
int SolverOptions::getNumberOfSweeps() const noexcept
{
    return static_cast<int> (pImpl->mSweeps);
}

/// Sets the convergence tolerance
void SolverOptions::setConvergenceTolerance(const double tolerance) noexcept
{
    pImpl->mConvergenceTolerance = tolerance;
}

/// Gets the convergence tolerance.
double SolverOptions::getConvergenceTolerance() const noexcept
{
    return pImpl->mConvergenceTolerance;
}

/// Sets the spherical solver radius in grid points
void SolverOptions::setSphericalSolverRadius(int radius) noexcept
{
    pImpl->mSphericalRadius = radius;
}

/// Gets the spherical solver radius in grid points.
int SolverOptions::getSphericalSolverRadius() const noexcept
{
    return pImpl->mSphericalRadius;
}

/// Sets the verbosity
void SolverOptions::setVerbosity(const EikonalXX::Verbosity verbosity) noexcept
{
    pImpl->mVerbosity = verbosity;
}

/// Gets the verbosity
EikonalXX::Verbosity SolverOptions::getVerbosity() const noexcept
{
    return pImpl->mVerbosity;
}

/// Sets the solver algorithm
void SolverOptions::setAlgorithm(const SolverAlgorithm algorithm) noexcept
{
    pImpl->mAlgorithm = algorithm;
}

/// Gets the solver algorithm
SolverAlgorithm SolverOptions::getAlgorithm() const noexcept
{
    return pImpl->mAlgorithm;
}

/// Resets the class and restores the defaults
void SolverOptions::clear() noexcept
{
    pImpl->mConvergenceTolerance = DEFFAULT_CONVERGENCE_TOLERANCE;
    pImpl->mSphericalRadius = DEFAULT_SPHERICAL_RADIUS;
    pImpl->mAlgorithm = SolverAlgorithm::LEVEL_SET_METHOD;
    pImpl->mVerbosity = Verbosity::ERROR;
    pImpl->mSweeps = DEFAULT_NUMBER_OF_SWEEPS;
}

