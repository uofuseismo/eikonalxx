#ifndef EIKONALXX_GEOMETRY_2D_HPP
#define EIKONALXX_GEOMETRY_2D_HPP
#include <memory>
namespace EikonalXX
{
/// @class Geometry2D "geometry2d.hpp" "eikonalxx/geometry2d.hpp"
/// @brief Defines the 2D physical model geometry.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Geometry2D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Geometry2D();
    /// @brief Copy constructor.
    /// @param[in] geometry  The geometry class from which to initialize this
    ///                      class.
    Geometry2D(const Geometry2D &geometry);
    /// @brief Move constructor.
    /// @param[in,out] geometry  The geometry class from which to initialize
    ///                          this class.  On exit, geometry's behavior
    ///                          is undefined.
    Geometry2D(Geometry2D &&geometry) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] geometry   The geometry class to copy to this.
    /// @result A deep copy of the geometry class.
    Geometry2D& operator=(const Geometry2D &geometry);
    /// @brief Move assignment operator.
    /// @param[in,out] geometry   The geometry class whose memory will be moved
    ///                           to this.  On exit, geometry's behavior is
    ///                           undefined.
    /// @result The memory from geometry moved to this.
    Geometry2D& operator=(Geometry2D &&geometry) noexcept; 
    /// @}

    /// @name Properties in x direction
    /// @{
    /// @brief Sets the number of grid points in x.
    /// @param[in] nx  The number of grid points in x.  This must be at least 3.
    /// @throws std::invalid_argument if nx is too small.
    /// @note This is required by the solver.
    void setNumberOfGridPointsInX(int nx);
    /// @result The number of grid points in x.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInX() is false.
    [[nodiscard]] int getNumberOfGridPointsInX() const;
    /// @result The number of cells in x.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInX() is false.
    [[nodiscard]] int getNumberOfCellsInX() const;
    /// @result True indicates that the number of grid points in x is set.
    [[nodiscard]] bool haveNumberOfGridPointsInX() const noexcept;

    /// @brief Sets the grid spacing in x.  This is required for the solver.
    /// @param[in] dx   The grid spacing in x in meters.
    /// @throws std::invalid_argument if dx is not positive.
    /// @note This is required by the solver. 
    void setGridSpacingInX(double dx);
    /// @result The grid spacing in x in meters.
    /// @throws std::runtime_error if \c haveGridSpacingInX() is false.
    [[nodiscard]] double getGridSpacingInX() const;
    /// @result True indicates that the grid spacing in x was set. 
    [[nodiscard]] bool haveGridSpacingInX() const noexcept;

    /// @brief Sets the model's x origin.
    /// @param[in] x0   The model's x origin in meters.
    /// @note By default this is 0.
    void setOriginInX(double x0) noexcept;
    /// @result The model's x origin in meters.
    double getOriginInX() const noexcept; 
    /// @}

    /// @name Properties in z direction
    /// @{
    /// @brief Sets the number of grid points in z.
    /// @param[in] nx  The number of grid points in z.  This must be at least 3.
    /// @throws std::invalid_argument if nz is too small.
    /// @note This is required by the solver.
    void setNumberOfGridPointsInZ(int nz);
    /// @result The number of grid points in z.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfGridPointsInZ() const;
    /// @result The number of cells in x.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfCellsInZ() const;
    /// @result True indicates that the number of grid points in z is set.
    [[nodiscard]] bool haveNumberOfGridPointsInZ() const noexcept;
 
    /// @brief Sets the grid spacing in z.
    /// @param[in] dz   The grid spacing in z in meters.
    /// @throws std::invalid_argument if dz is not positive.
    /// @note This is required by the solver.
    void setGridSpacingInZ(double dz);
    /// @result The grid spacing in z in meters.
    /// @throws std::runtime_error if \c haveGridSpacingInZ() is false.
    [[nodiscard]] double getGridSpacingInZ() const;
    /// @result True indicates that the grid spacing in z was set.
    [[nodiscard]] bool haveGridSpacingInZ() const noexcept;

    /// @brief Sets the model's z origin.
    /// @param[in] z0   The model's z origin in meters.
    /// @note By default this is 0.
    void setOriginInZ(double z0) noexcept;
    /// @result The model's z origin in meters.
    double getOriginInZ() const noexcept;
    /// @}

    /// @result The number of grid points in the geometry.
    /// @throws std::runtime_error if \c haveNumberOfGridPointsInX() is false
    ///         or \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfGridPoints() const;
    /// @result The number of cells in the geometry.
    /// @throws std::runtime_error if \c haveNumberOfGridPointsInX() is false
    ///         or \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfCells() const;

    /// @name Destructors
    /// @{
    /// @brief Destructor.
    ~Geometry2D();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @} 
private:
    class Geometry2DImpl;
    std::unique_ptr<Geometry2DImpl> pImpl;
};
/// @result True if the left-hand side geometry is equivalent to the
///         right-hand side geometry.
bool operator==(const Geometry2D &lhs, const Geometry2D &rhs);
/// @result True if the left-hand side geometry is not equivalent to
///         the right-hand side geometry.
bool operator!=(const Geometry2D &lhs, const Geometry2D &rhs);
}
#endif
