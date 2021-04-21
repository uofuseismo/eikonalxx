#ifndef PYEIKONALXX_SOLVEROPTIONS_HPP
#define PYEIKONALXX_SOLVEROPTIONS_HPP
#include <memory>
#include <pybind11/pybind11.h>
#include <eikonalxx/enums.hpp>
namespace EikonalXX
{
class SolverOptions;
}
namespace PEikonalXX
{
class SolverOptions
{
public:
    /// @brief Constructor.
    SolverOptions();
    /// @brief Copy c'tor.
    SolverOptions(const SolverOptions &options);
    SolverOptions(const EikonalXX::SolverOptions &options);
    /// @brief Copy c'tor.
    SolverOptions(SolverOptions &&options) noexcept;
    /// @brief Copy assignment.
    SolverOptions& operator=(const SolverOptions &options);
    SolverOptions& operator=(const EikonalXX::SolverOptions &options);
    /// @brief Move assignment.
    SolverOptions& operator=(SolverOptions &&options) noexcept;
    /// @brief Destructor
    ~SolverOptions();
    /// @brief Reset teh class
    void clear() noexcept;

    void setNumberOfSweeps(uint16_t nSweeps) noexcept;
    /// @result Gets the number of Gauss-Seidel sweeps.
    [[nodiscard]] int getNumberOfSweeps() const noexcept;

    /// @brief Sets the convergence tolerance.
    void setConvergenceTolerance(double tolerance) noexcept; 
    /// @result The convergence tolerance in seconds.
    [[nodiscard]] double getConvergenceTolerance() const noexcept;

    /// @brief Sets the spherical to cartesian solver transition distance.
    void setSphericalSolverRadius(int epsilon) noexcept;
    /// @result The radius in grid points from the source where the
    ///         spherical finite-difference stencils are employed.
    [[nodiscard]] int getSphericalSolverRadius() const noexcept;

    /// @brief Sets the solver verbosity.
    /// @param[in] verbosity  The solver verbosity level.
    void setVerbosity(EikonalXX::Verbosity verbosity) noexcept;
    /// @result The solver verbosity.
    [[nodiscard]] EikonalXX::Verbosity getVerbosity() const noexcept;

    /// @brief Defines the solver algorithm. 
    /// @param[in] algorithm  The solver algorithm.
    void setAlgorithm(EikonalXX::SolverAlgorithm algorithm) noexcept;
    /// @result The solver algorithm.
    [[nodiscard]] EikonalXX::SolverAlgorithm getAlgorithm() const noexcept;
private:
    std::unique_ptr<EikonalXX::SolverOptions> pImpl;
};
void initializeSolverOptions(pybind11::module &m);
}
#endif
