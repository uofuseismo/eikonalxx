#ifndef EIKONALXX_RAY_TRACER_2D_HPP
#define EIKONALXX_RAY_TRACER_2D_HPP
#include <vector>
#include <memory>
namespace EikonalXX
{
 namespace AbstractBaseClass
 {
  template<class T> class ISolver2D;
 }
 class Geometry2D;
 class Station2D;
 namespace Ray
 {
  class Path2D;
 }
}
namespace EikonalXX::Ray
{
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

    /// @name Step 3: Perform Ray Tracing
    /// @{

    /// @brief Traces the rays through the travel time field.
    /// @throws std::invalid_argument if \c solver.haveGradientTravelTimeField()
    ///         is false.
    /// @throws std::runtime_error if \c haveStations() is false.
    template<typename T>
    void trace(const EikonalXX::AbstractBaseClass::ISolver2D<T> &solver);
    /// @}

    /// @name Step 4: Fetch Results
    /// @{

    /// @result True indicates the ray paths were trace.
    [[nodiscard]] bool haveRayPaths() const noexcept;
    /// @result The source-to-receiver ray paths for each station.
    /// @throws std::runtime_error if \c haveRayPaths() is false.
    [[nodiscard]] std::vector<Path2D> getRayPaths() const;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] tracer  The tracer to copy to this.
    /// @result A deep copy of the input tracer.
    GradientTracer2D& operator=(const GradientTracer2D &tracer);
    /// @brief Move assignment.
    /// @param[in] tracer  The tracer whose memory will be moved to this.
    ///                    On exit, tracer's behavior is undefined.
    /// @result the memory from tracer moved to this.
    GradientTracer2D& operator=(GradientTracer2D &&tracer) noexcept;
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
