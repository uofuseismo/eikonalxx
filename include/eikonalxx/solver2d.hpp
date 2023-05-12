#ifndef EIKONALXX_SOLVER_2D_HPP
#define EIKONALXX_SOLVER_2D_HPP
#include <memory>
#include <vector>
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
namespace EikonalXX
{
// Forward declarations
class Geometry2D;
class Source2D;
template<class T> class Model2D;
class SolverOptions;
/// @class Solver2D "solver2d.hpp" "eikonalx/solver2d.hpp"
/// @brief Solves the eikonal equation in a 2D slice.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class Solver2D : public EikonalXX::AbstractBaseClass::ISolver2D<T>
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Solver2D();
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
                    const Geometry2D &geometry);
    /// @result The solver options.
    /// @throws std::runtime_erorr if \c isInitialized() is false. 
    [[nodiscard]] SolverOptions getOptions() const;
    /// @result The model's geometry.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] Geometry2D getGeometry() const;
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
    void setVelocityModel(const Model2D<T> &velocityModel);
    /// @result The velocity model.
    /// @throws std::runtime_error if the \c haveVelocityModel() is false.
    [[nodiscard]] Model2D<T> getVelocityModel() const;
    /// @result True indicates that the velocity model was set.
    [[nodiscard]] bool haveVelocityModel() const noexcept;
    /// @}

    /// @name Step 3: Source
    /// @{

    /// @brief Sets the source location.
    /// @param[in] source   A class defining the source.
    /// @throws std::invalid_argument if the source location in x and z
    ///         is not set.
    /// @throws std::runtime_error if the class is not initialized.
    void setSource(const Source2D &source);
    /// @brief Sets the source location.
    /// @param[in] location   The source location.  location.first is the x
    ///                       position in the model and location.second is
    ///                       the z location in the model.
    /// @throws std::runtime_error if the class is not initialized.
    /// @throws std::invalid_argument if the source location is not in the
    ///         model.
    /// @sa \c isInitialized()
    void setSource(const std::pair<double, double> &location);
    /// @result The source information.
    /// @throws std::runtime_error if \c haveSource() is false.
    [[nodiscard]] Source2D getSource() const;
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

    /// @name Step 5: Results
    /// @{

    /// @result The travel times from the source to all nodes in the model in
    ///         in seconds.  This uses the natural ordering.
    /// @note This has dimension getGeometry.getNumberOfGridPoints().
    /// @sa \c Ordering2D, \c haveTravelTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeField() const override;
    /// @result A pointer to the travel time field at all nodes in the model
    ///         in seconds.   This uses the natural ordering and
    ///         has dimension [getGeometry.getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c haveTravelTimeField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T* getTravelTimeFieldPointer() const override;
    /// @result True indicates that \c solve() has been called and the travel
    ///         time field is available.
    [[nodiscard]] bool haveTravelTimeField() const noexcept override;

    /// @result The gradient of the travel time field in x in seconds/meter.
    ///         This uses the natural ordering.
    /// @note This has dimension getGeometry.getNumberOfGridPoints().
    /// @sa \c Ordering2D, \c haveTravelTimeGradientField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeGradientFieldInX() const override;
    /// @result A reference to the gradient of the travel time field in x
    ///         in seconds/meter.  This uses the natural ordering and
    ///         has dimension [getGeometry.getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c haveTravelTimeGradientField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T *getTravelTimeGradientFieldInXPointer() const override;

    /// @result The gradient of the travel time field in z in seconds/meter.
    ///         This uses the natural ordering.
    /// @note This has dimension getGeometry.getNumberOfGridPoints().
    /// @sa \c Ordering2D, \c haveTravelGradientTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeGradientFieldInZ() const override;
    /// @result A reference to the gradient of the travel time field in z
    ///         in seconds/meter.  This uses the natural ordering and
    ///         has dimension [getGeometry.getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c haveTravelTimeGradientField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T *getTravelTimeGradientFieldInZPointer() const override;

    /// @result True indicates the travel time field's gradient was computed.
    [[nodiscard]] bool haveTravelTimeGradientField() const noexcept override;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Solver2D() override;
    /// @}

    // Remove some functionality.
    Solver2D(const Solver2D &solver) = delete;
    Solver2D(Solver2D &&solver) noexcept = delete;
    Solver2D& operator=(const Solver2D &solver) = delete;
    Solver2D& operator=(Solver2D &&solver) noexcept = delete;
private:
    class Solver2DImpl;
    std::unique_ptr<Solver2DImpl> pImpl;
};
}
#endif
