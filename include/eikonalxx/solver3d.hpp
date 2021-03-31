#ifndef EIKONALXX_SOLVER3D_HPP
#define EIKONALXX_SOLVER3D_HPP
#include <memory>
#include <vector>
namespace EikonalXX
{
// Forward declarations
class Geometry3D;
class Source3D;
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

    /// @name Step 2: Velocity Model
    /// @{
    /// @brief Sets the velocity model.
    /// @param[in] velocityModel  The velocity model.
    /// @throws std::runtime_error if the class is not initialized.
    /// @throws std::invalid_argument if the velocity model was not
    ///         set on velocityModel or its dimensions are inconsistent with
    ///         the geometry set during initialization.
    /// @sa \c isInitialized()
    void setVelocityModel(const Model3D<T> &velocityModel);
    /// @result The velocity model.
    /// @throws std::runtime_error if the \c haveVelocityModel() is false.
    [[nodiscard]] Model3D<T> getVelocityModel() const;
    /// @result True indicates that the velocity model was set.
    [[nodiscard]] bool haveVelocityModel() const noexcept;
    /// @}

    /// @name Step 3: Source
    /// @{
    /// @brief Sets the source location.
    /// @param[in] source   A class defining the source.
    /// @throws std::invalid_argument if the source location in x, y, and z
    ///         is not set.
    /// @throws std::runtime_error if the class is not initialized.
    void setSource(const Source3D &source);
    /// @brief Sets the source location.
    /// @param[in] location   The source location.  std::get<0> (location)
    ///                       is the x position in the model,
    ///                       std::get<1> (location) is the y position in the
    ///                       model, and std::get<2> (location) is the z
    ///                       location in the model.
    /// @throws std::runtime_error if the class is not initialized.
    /// @throws std::invalid_argument if the source location is not in the
    ///         model.
    /// @sa \c isInitialized()
    void setSource(const std::tuple<double, double, double> &location);
    /// @result The source information.
    /// @throws std::runtime_error if \c haveSource() is false.
    [[nodiscard]] Source3D getSource() const;
    /// @result True indicates that the source was set.
    [[nodiscard]] bool haveSource() const noexcept;
    // @}

    /// @name Step 4: Solve
    /// @{
    /// @brief Solves the eikonal equation for the given source/velocity model.
    /// @throws std::runtime_error if the source or velocity model is not set.
    /// @sa \c isInitialized(), \c haveVelocityModel(), \c haveSource()
    void solve(); 
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
