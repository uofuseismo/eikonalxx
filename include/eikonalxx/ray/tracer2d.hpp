#ifndef EIKONALXX_RAY_TRACER_2D_HPP
#define EIKONALXX_RAY_TRACER_2D_HPP
#include <vector>
#include <memory>
namespace EikonalXX
{
 template<class T> class Solver2D;
 class Station2D;
 namespace Ray
 {
  class RayPath2D;
 }
}
namespace EikonalXX::Ray
{
template<class T>
/// @class GradientTracer2D "gradientTracer2d.hpp" "eikonalxx/ray/gradientTrace2d.hpp"
/// @brief This class computes ray paths by marching up the travel time field's
///        gradient from the receiver to the source.   
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class GradientTracer2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    GradientTracer2D();
    /// @brief Copy constructor.
    /// @param[in] tracer  The tracer class from which to initialize this class.
    GradientTracer2D(const GradientTracer2D &tracer);
    /// @brief Move constructor.
    /// @param[in,out] tracer  The tracer class from which to initialize this
    ///                        class.  On exit, tracer's behavior is undefined.
    GradientTracer2D(GradientTracer2D &&tracer) noexcept;
    /// @}
 
    /// @name Step 1: Initialization
    /// @{
 
    /// @brief Sets the geometry and options.
    void initialize(const EikonalXX::Geometry2D &geometry); 
    /// @result True indicates the class is initialized. 
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @}
 
    /// @name Step 2: Set Stations
    /// @{

    /// @brief Sets the stations at which to compute the rays.
    /// @throws std::invalid_argument if the any stations are not in the
    ///         computational grid.
    void setStations(const std::vector<EikonalXX::Station2D> &stations);
    /// @result True indicates that stations were set.
    [[nodiscard]] bool haveStations() const noexcept;
    /// @}

    /// @name Step 3: Set Travel Time Fields
    /// @{

    /// @brief Sets the travel time field.
    /// @param[in] solver  The solver containing the travel time field.
    /// @throws std::invalid_argument if there is no source position on the
    ///         solver, the travel times were not yet computed, or the
    ///         geometries are inconsistent.
    /// @throws std::runtime_error if \c isInitialized() is false.
    void setTravelTimeField(const EikonalXX::Solver2D<T> &solver);
    /// @result True indicates the travel time field was set.
    [[nodiscard]] bool haveGradientTravelTimeFields() const noexcept;
    /// @}

    /// @name Step 4: Perform Ray Tracing
    /// @{

    /// @brief Traces the rays through the travel time field.
    /// @throws std::runtime_error if \c haveGradientTravelTimeFields()
    ///         is false.
    void trace();
    /// @}

    /// @name Step 5: Fetch Results
    /// @{

    /// @result True indicates the ray paths were trace.
    [[nodiscard]] bool haveRayPaths() const noexcept;
    /// @result The source-to-receiver ray paths for each station.
    /// @throws std::runtime_error if \c haveRayPaths() is false.
    [[nodiscard]] std::vector<RayPath2D> getRayPaths() const;
    /// @}

    /// @name Operators
    /// @{
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~GradientTracer2D();
    /// @}
private:
    class GradientTracer2DImpl;
    std::unique_ptr<GradientTracer2DImpl> pImpl;
};
}
#endif
