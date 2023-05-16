#ifndef EIKONALXX_ANALYTIC_HOMOGENEOUS_3D_HPP
#define EIKONALXX_ANALYTIC_HOMOGENEOUS_3D_HPP
#include <memory>
#include <vector>
#include "eikonalxx/abstractBaseClass/solver3d.hpp"
namespace EikonalXX
{
/// Forward declarations
class Geometry3D;
class Source3D;
namespace Analytic
{
/// @class Homogeneous3D "homogeneous3d.hpp" "eikonalxx/analytic/homogeneous3d.hpp"
/// @brief Solves the eikonal equation in a constant velocity model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
/// @ingroup Solver_Analytic_chapter
template<class T>
class Homogeneous3D : public EikonalXX::AbstractBaseClass::ISolver3D<T>
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Homogeneous3D(); 
    /// @brief Copy constructor.
    Homogeneous3D(const Homogeneous3D &solver);
    /// @brief Move constructor.
    Homogeneous3D(Homogeneous3D &&solver) noexcept;
    /// @}

    /// @name Step 1: Initialization
    /// @{

    /// @brief Initializes the class.
    void initialize(const EikonalXX::Geometry3D &geometry);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @result The model geometry. 
    /// @throws std::runtime_error if \c isInitialized() is false. 
    [[nodiscard]] Geometry3D getGeometry() const;
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
    /// @}

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
    /// @sa \c Ordering3D, \c haveTravelTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeField() const override;
    /// @result A pointer to the travel time field at all nodes in the model
    ///         in seconds.   This uses the natural ordering and
    ///         has dimension [getGeometry.getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c haveTravelTimeField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering3D
    [[nodiscard]] const T* getTravelTimeFieldPointer() const override;
    /// @result True indicates that \c solve() has been called and the travel
    ///         time field is available.
    [[nodiscard]] bool haveTravelTimeField() const noexcept override;

    /// @result The gradient of the travel time field in x, y, and z - all in
    ///         seconds/meter.  This uses the natural ordering.  This has
    ///         dimension [getGeometry.getNumberOfGridPoints() x 3] and is
    ///         stored in row major format.
    /// @sa \c Ordering2D, \c haveTravelTimeField(), \c getGeometry()
    [[nodiscard]] std::vector<T> getTravelTimeGradientField() const override;
    /// @result A pointer to the gradient of the travel time field in x, y, 
    ///         and z - all in seconds/meter.  This uses the natural ordering.
    ///         This has dimension [getGeometry.getNumberOfGridPoints() x 3]
    ///         and is stored in row major format.
    /// @throws std::runtime_error if \c haveTravelTimeGradientField() is false.
    /// @sa \c haveTravelTimeField(), \c getGeometry(), \c Ordering2D
    [[nodiscard]] const T *getTravelTimeGradientFieldPointer() const override;

    /// @result True indicates the travel time field's gradient was computed.
    [[nodiscard]] bool haveTravelTimeGradientField() const noexcept override;

    /// @brief Writes the travel time field to VTK.
    /// @param[in] fileName  The name of the VTK file.
    /// @param[in] title     The dataset's title.
    /// @throws std::runtime_error if \c haveTravelTimeField() is false.
    /// @throws std::invalid_argument if there is an error while opening the
    ///         output file.
    void writeVTK(const std::string &fileName,
                  const std::string &title = "homogeneous_analytic_traveltime_field") const;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] solver  The 3D homogeneous analytic solver to copy to this.
    /// @result A deep copy of the input solver.
    Homogeneous3D& operator=(const Homogeneous3D &solver);
    /// @brief Move assignment.
    /// @param[in,out] solver  The 3D homogeneous analytic solver whose memory
    ///                        will be moved to this.  On exit, solver's
    ///                        behavior is undefined.
    /// @result The memory from the input solver moved to this.
    Homogeneous3D& operator=(Homogeneous3D &&solver) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Releases all memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Homogeneous3D() override;
    /// @}
private:
    class Homogeneous3DImpl;
    std::unique_ptr<Homogeneous3DImpl> pImpl;
};
}
}
#endif
