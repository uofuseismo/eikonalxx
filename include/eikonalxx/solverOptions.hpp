#ifndef EIKONALXX_SOLVER_OPTIONS_HPP
#define EIKONALXX_SOLVER_OPTIONS_HPP
#include <ostream>
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
/// @class SolverOptions "solverOptions.hpp" "eikonalxx/solverOptions.hpp"
/// @brief Defines the fast-sweeping method eikonal solver options.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class SolverOptions
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    SolverOptions();
    /// @brief Copy constructor.
    /// @param[in] options   The options class from which to initialize
    ///                       this class.
    SolverOptions(const SolverOptions &options);
    /// @brief Move constructor.
    /// @param[in,out] options  The options class from which to initialize
    ///                         this class.  On exit, options's behavior is
    ///                         undefined. 
    SolverOptions(SolverOptions &&options) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] options   The options class to copy to this.
    /// @result A deep copy of the options class.
    SolverOptions& operator=(const SolverOptions &options);
    /// @brief Move assignment operator.
    /// @param[in,out] options  The options class whose memory will be moved
    ///                         to this.  On exit, options's behavior is
    ///                         undefined.
    /// @result The memory from options moved to this.
    SolverOptions& operator=(SolverOptions &&options) noexcept;
    /// @} 

    /// @name Properties
    /// @{

    /// @brief Sets the number of (Gauss-Seidel) sweeps.
    /// @param[in] nSweeps  The number of Gauss-Seidel sweeps.
    /// @note This does not include the initialization sweep which will always
    ///       be performed.  Additionally, the solver may convergence
    ///       before the number of specified sweeps is performed.
    /// @sa \c setConvergenceTolerance()
    void setNumberOfSweeps(uint16_t nSweeps) noexcept;
    /// @result Gets the number of Gauss-Seidel sweeps.
    [[nodiscard]] int getNumberOfSweeps() const noexcept;

    /// @brief Sets the convergence tolerance of the underlying iterative
    ///        Gauss-Seidel method.
    /// @param[in] tolerance  The convergence tolerance in seconds.
    ///                       The iterative solver with terminate early if
    ///                       if the max change in travel time at a node
    ///                       is less than this value.  Setting tolerance
    ///                       to zero or a negative number effectively
    ///                       disables it.
    /// @sa \c setNumberOfSweeps()
    void setConvergenceTolerance(double tolerance) noexcept; 
    /// @result The convergence tolerance in seconds.
    [[nodiscard]] double getConvergenceTolerance() const noexcept;

    /// @brief Near the source the wavefront has high curvature and the 
    ///        Cartesian finite-difference operator is inaccurate.  
    ///        To mitigate this a solver using the factored eikonal
    ///        eikonal equation is used to update the traveltimes
    ///        near the source.
    /// @param[in] epsilon  The number of grid points the update grid point
    ///                     must be from the source grid point.  This 
    ///                     effectively draws a box around the source grid
    ///                     point.  Setting this to less than 0 means that
    ///                     the Cartesian solver will be used for all updates.
    void setFactoredEikonalEquationSolverRadius(int epsilon) noexcept;
    /// @result The number of grid points in x and z in 2D,
    ///         or x, y, and z in 3D that the update grid point
    ///         must be from the source grid point to use the Cartesian solver.
    [[nodiscard]] int getFactoredEikonalEquationSolverRadius() const noexcept;

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
    ///@}

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~SolverOptions();
    /// @brief Resets the class and restores the defaults.
    void clear() noexcept;
    /// @}
private:
    class SolverOptionsImpl;
    std::unique_ptr<SolverOptionsImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const SolverOptions &options);
}
#endif
