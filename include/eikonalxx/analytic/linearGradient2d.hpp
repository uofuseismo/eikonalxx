#ifndef EIKONALXX_ANALYTIC_LINEAR_GRADIENT_2D_HPP
#define EIKONALXX_ANALYTIC_LINEAR_GRADIENT_2D_HPP
#include <memory>
#include <vector>
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
namespace EikonalXX
{
/// Forward declarations
class Geometry2D;
class Source2D;
namespace Analytic
{
/// @class LinearGradient2D "linearGradient2d.hpp" "eikonalxx/analytic/linearGradient2d.hpp"
/// @brief Solves the eikonal equation in a velocity model with a
///        linear gradient.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @note Further details can be found in Fomel's : Fast sweeping method for
///       the factored eikonal equation.
/// @ingroup Solver_Analytic_chapter
template<class T>
class LinearGradient2D : public EikonalXX::AbstractBaseClass::ISolver2D<T>
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    LinearGradient2D(); 
    /// @brief Copy constructor.
    LinearGradient2D(const LinearGradient2D &solver);
    /// @brief Move constructor.
    LinearGradient2D(LinearGradient2D &&solver) noexcept;
    /// @}

    /// @name Step 1: Initialization
    /// @{

    /// @brief Initializes the class.
    void initialize(const EikonalXX::Geometry2D &geometry);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The model geometry. 
    /// @throws std::runtime_error if \c isInitialized() is false. 
    [[nodiscard]] Geometry2D getGeometry() const;
    /// @}
    
    /// @name Step 2: Velocity Model
    /// @{

    /// @brief Sets the linear gradient velocity model.
    /// @param[in] velocity   velocity.first is the velocity in m/s at the
    ///                       top of the model and velocity.second is the
    ///                       velocity in m/s at the bottom of the model.
    /// @throws std::invalid_argument if velocity is not positive.
    void setVelocityModel(const std::pair<double, double> &velocity);
    /// @result True indicates that the velocity model was set.
    [[nodiscard]] bool haveVelocityModel() const noexcept;
    /// @result The slowness in s/m at a cell.
    /// @param[in] iCellX  The cell index in x.
    /// @param[in] iCellZ  The cell index in z.
    /// @throws std::runtime_error if \c haveVelocityModel() is false.
    /// @throws std::invalid_argument if iCellX or iCellZ is out of the geometry
    ///         bounds.
    [[nodiscard]] T getSlowness(const int iCellX, const int iCellZ) const override;
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
    [[nodiscard]] Source2D getSource() const override;
    /// @result True indicates that the source was set.
    [[nodiscard]] bool haveSource() const noexcept override;
    /// @}

    /// @name Step 4: Solve
    /// @{

    /// @brief Solves the eikonal equation for the given source/velocity model.
    /// @throws std::runtime_error if the source or velocity model is not set.
    /// @sa \c isInitialized(), \c haveVelocityModel(), \c haveSource()
    void solve();
    /// @brief Computes the gradient travel time field.
    /// @throws std::runtime_error if the source or velocity model is not set.
    /// @sa \c isInitialized(), \c haveVelocityModel(), \c haveSource()
    void computeTravelTimeGradientField();
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

    /// @result The gradient of the travel time field in x and z - both in
    ///         seconds/meter.  This uses the natural ordering.  This has
    ///         dimension [getGeometry.getNumberOfGridPoints() x 2] and is
    ///         stored in row major format.
    /// @sa \c Ordering2D, \c haveTravelTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeGradientField() const override;
    /// @result A pointer to the gradient of the travel time field in x and z -
    ///         both in seconds/meter.  This uses the natural ordering.  This
    ///         has dimension [getGeometry.getNumberOfGridPoints() x 2] and
    ///         is stored in row major format.
    /// @throws std::runtime_error if \c haveTravelTimeGradientField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T *getTravelTimeGradientFieldPointer() const override;
    /// @result True indicates the gradient of travel time field was computed.
    [[nodiscard]] bool haveTravelTimeGradientField() const noexcept override;

    /// @brief Writes the travel time field to VTK.
    /// @param[in] fileName  The name of the VTK file.
    /// @param[in] title     The dataset's title.
    /// @param[in] writeGradient  If true then also write the gradient of the
    ///                           travel time field.
    /// @throws std::runtime_error if \c haveTravelTimeField() is false.
    ///         Additionally, this throws if writeGradient is true and
    ///         \c haveGradientTravelTimefield() is false.
    /// @throws std::invalid_argument if there is an error while opening the
    ///         output file.
    void writeVTK(const std::string &fileName,
                  const std::string &title = "linearGradient_analytic_traveltime_field",
                  const bool writeGradient = false) const;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] solver  The 2D linear graident analytic solver to copy
    ///                    to this.
    /// @result A deep copy of the input solver.
    LinearGradient2D& operator=(const LinearGradient2D &solver);
    /// @brief Move assignment.
    /// @param[in,out] solver  The 2D linear gradient analytic solver whose
    ///                        memory will be moved to this.  On exit, solver's
    ///                        behavior is undefined.
    /// @result The memory from the input solver moved to this.
    LinearGradient2D& operator=(LinearGradient2D &&solver) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Releases all memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~LinearGradient2D() override;
    /// @}
private:
    class LinearGradient2DImpl;
    std::unique_ptr<LinearGradient2DImpl> pImpl;
};
}
}
#endif
