#ifndef EIKONALXX_ANALYTIC_HOMOGENEOUS_2D_HPP
#define EIKONALXX_ANALYTIC_HOMOGENEOUS_2D_HPP
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
/// @class Homogeneous2D "homogeneous2d.hpp" "eikonalxx/analytic/homogeneous2d.hpp"
/// @brief Solves the eikonal equation in a constant velocity model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T>
class Homogeneous2D : public EikonalXX::AbstractBaseClass::ISolver2D<T>
{
public:
    /// @name Constructor
    /// @{

    /// @brief Constructor.
    Homogeneous2D(); 
    /// @brief Copy constructor.
    Homogeneous2D(const Homogeneous2D &solver);
    /// @brief Move constructor.
    Homogeneous2D(Homogeneous2D &&solver) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] solver  The 2D homogeneous analytic solver to copy to this.
    /// @result A deep copy of the input solver.
    Homogeneous2D& operator=(const Homogeneous2D &solver);
    /// @brief Move assignment.
    /// @param[in,out] solver  The 2D homogeneous analytic solver whose memory
    ///                        will be moved to this.  On exit, solver's
    ///                        behavior is undefined.
    /// @result The memory from the input solver moved to this.
    Homogeneous2D& operator=(Homogeneous2D &&solver) noexcept;
    /// @}

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

    /// @brief Sets the constant velocity model.
    /// @param[in] velocity   The velocity in m/s.
    /// @throws std::invalid_argument if velocity is not positive.
    void setVelocityModel(double velocity);
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
    /// @}

    /// @name Step 4: Solve
    /// @{

    /// @brief Solves the eikonal equation for the given source/velocity model.
    /// @throws std::runtime_error if the source or velocity model is not set.
    /// @sa \c isInitialized(), \c haveVelocityModel(), \c haveSource()
    void solve();
    /// @brief Computes the gradient travel time field.
    /// @throws std::runtime_error if the source or velocity model is not set.
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

    /// @result The gradient of the travel time field in x in seconds/meter.
    ///         This uses the natural ordering.
    /// @note This has dimension getGeometry.getNumberOfGridPoints().
    /// @sa \c Ordering2D, \c haveTravelTimeField(), \c getGeometry()
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
    /// @sa \c Ordering2D, \c haveTravelTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeGradientFieldInZ() const override;
    /// @result A reference to the gradient of the travel time field in z
    ///         in seconds/meter.  This uses the natural ordering and
    ///         has dimension [getGeometry.getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c haveTravelTimeGradientField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T *getTravelTimeGradientFieldInZPointer() const override;

    /// @result True indicates the travel time field's gradient was computed.
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
                  const std::string &title = "homogeneous_analytic_traveltime_field",
                  bool writeGradient = false) const;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Releases all memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Homogeneous2D() override;
    /// @}
private:
    class Homogeneous2DImpl;
    std::unique_ptr<Homogeneous2DImpl> pImpl;
};
}
}
#endif
