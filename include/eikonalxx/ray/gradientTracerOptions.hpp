#ifndef EIKONALXX_RAY_GRADIENT_TRACER_OPTIONS_HPP
#define EIKONALXX_RAY_GRADIENT_TRACER_OPTIONS_HPP
#include <vector>
#include <memory>
namespace EikonalXX::Ray
{
/// @class GradientTracerOptions "gradientTracerOptions.hpp" "eikonalxx/ray/gradientTracerOptions.hpp"
/// @brief Defines the gradient ray tracer options.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class GradientTracerOptions
{
public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    GradientTracerOptions();
    /// @brief Copy constructor.
    /// @param[in] options  The options class from which to initialize
    ///                     this class.
    GradientTracerOptions(const GradientTracerOptions &options);
    /// @brief Move constructor.
    /// @param[in,out] options  The options class from which to initialize this
    ///                         class.  On exit, options's behavior is
    ///                         undefined.
    GradientTracerOptions(GradientTracerOptions &&options) noexcept;
    /// @}

    /// @brief Sets the step length for a given distance from the source.
    ///        Effectively, when the ray's current cell in x, y, or z
    ///        is less than the given radius then the corresponding step
    ///        length in the gradient marching will be employed.
    ///        For completeness, every step along the gradient looks like:
    ///        \f[
    ///           \textbf{x}_1 = \textbf{x}_0
    ///                        + \alpha \min\{\Delta x, \Delta y, \Delta z\}
    ///                          \frac{\textbf{g}}{|\textbf{g}|}
    ///        \f]
    ///        where the scale factors are \f$ \alpha \f$ and \f$ \Delta x \f$,
    ///        \f$ \Delta y \f$, and \f$ \Delta z \f$ is the grid spacing in
    ///        each direction.
    /// @param[in] radiusScaleFactor  Each element of this vector defines a
    ///                               search radius (i.e., the distance of the
    ///                               ray's current cell to the source cell)
    ///                               and corresponding gradient scale factor.
    ///                               The scale factor must be positive.
    /// @throws std::invalid_argument if any distances are duplicated, any scale
    ///         factor is not positive, or any radius is negative.
    /// @note While you may be tempted to make one very small stepsize you will
    ///       end up with a lot of intermediate ray segments.
    void setRadiusScaleFactor(const std::vector<std::pair<int, double>> &radiusScaleFactor);
    /// @result The gradient scale factors for given distances sorted in 
    ///         ascending order of search radius.
    /// @note By default we take smaller steps as we get closer to the source
    ///       so we don't blow past it.
    [[nodiscard]] std::vector<std::pair<int, double>> getRadiusScaleFactor() const noexcept;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] options  The options class to copy to this.
    /// @result A deep copy of the options.
    GradientTracerOptions& operator=(const GradientTracerOptions &options);
    /// @brief Move assignment.
    /// @param[in,out] options  The options whose memory will be moved to this.
    ///                         On exit, options's behavior is undefined.
    /// @result The memory from options moved to this.
    GradientTracerOptions& operator=(GradientTracerOptions &&options) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and restores defaults.
    void clear() noexcept;
    /// @brief Destructor.
    ~GradientTracerOptions();
    /// @}
private:
    class GradientTracerOptionsImpl;
    std::unique_ptr<GradientTracerOptionsImpl> pImpl;
};
}
#endif
