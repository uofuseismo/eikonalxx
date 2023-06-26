#ifndef EIKONALXX_RAY_LAYER_SOLVER_HPP
#define EIKONALXX_RAY_LAYER_SOLVER_HPP
#include <memory>
namespace EikonalXX::Ray
{
class Segment2D;
}
namespace EikonalXX::Ray
{
/// @brief A simple utility class for computing rays for a first arrival
///        in 1D, isotropic, layered media.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class LayerSolver
{
public:
    /// @brief Constructor.
    LayerSolver();
 
    /// @brief Sets a layer cake model of the form:
    ///        Interface 0 ---------------------------- (Free Surface)
    ///                    |        Velocity 0        |
    ///        Interface 1 ----------------------------
    ///                    |        Velocity 1        |
    ///                    ----------------------------
    ///                                 .
    ///                                 .
    ///                                 .
    ///        Interface n --------------------------- (Underlying half space)
    ///                    |        Velocity n       |
    ///                    ---------------------------
    ///
    /// @param[in] interfaces     The interfaces for each layer.
    /// @param[in] velocityModel  The velocities in each layer.
    /// @note Velocity inversions are not yet considered.
    void setVelocityModel(const std::vector<double> &interfaces,
                          const std::vector<double> &velocityModel);
    /// @result True indicates the velocity model was set.
    [[nodiscard]] bool haveVelocityModel() const noexcept;

    /// @brief Sets the source depth.
    /// @param[in] depth  The source depth in meters.
    /// @throws std::invalid_argument if the source is not in the layer cake
    ///         model.
    /// @throws std::runtime_error if \c haveVelocityModel() is false.
    void setSourceDepth(const double depth);
    /// @result The source depth in meters.
    [[nodiscard]] double getSourceDepth() const;
    /// @result True indicates the source depth was set.
    [[nodiscard]] bool haveSourceDepth() const noexcept;

    /// @brief Sets the station offset, in meters, with the station at the
    ///        free surface.
    void setStationOffset(double offset);
    /// @brief Sets the station offset, in meters, and the station depth
    ///        in meters.
    void setStationOffsetAndDepth(double offset, double depth);
    [[nodiscard]] bool haveStationOffsetAndDepth() const noexcept;

    void solve();

    /// @brief Reset the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~LayerSolver();
private:
    class LayerSolverImpl;
    std::unique_ptr<LayerSolverImpl> pImpl;
};
}
#endif
