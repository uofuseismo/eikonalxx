#ifndef EIKONALXX_SOLVER3D_HPP
#define EIKONALXX_SOLVER3D_HPP
#include <memory>
#include <vector>
namespace EikonalXX
{
// Forward declarations
class Geometry3D;
template<class T> class Model3D;
class SolverOptions;
/// @class Solver3D "solver3d.hpp" "eikonalx/solver3d.hpp"
/// @brief Solves the eikonal equation in a 3D medium.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class Solver3D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Solver3D();
    /// @}

    /// @name Step 1: Initialization
    /// @{
    /// @brief Initializes the class.
    /// @param[in] options   Contains the solver options.
    /// @param[in] geometry  Defines the model geometry. 
    /// @throws std::invalid_argument if the grid spacing or number of
    ///         grid points in all directions is not set.
    /// @note There is a substantial amount of allocation and preprocessing
    ///       that occurs in this function so it is best to call it sparingly.
    void initialize(const SolverOptions &options,
                    const Geometry3D &geometry);
    /// @result The solver options.
    /// @throws std::runtime_erorr if \c isInitialized() is false. 
    [[nodiscard]] SolverOptions getOptions() const;
    /// @result The model's geometry.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] Geometry3D getGeometry() const;
    /// @result True indicates that the class is initialized. 
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @}

    /// @brief Destructor.
    ~Solver3D();
    /// @brief Resets the class and restores defaults.
    void clear() noexcept;
private:
    class Solver3DImpl;
    std::unique_ptr<Solver3DImpl> pImpl;
};
}
#endif
