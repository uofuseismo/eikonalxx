#ifndef EIKONALXX_GEOMETRY3D_HPP
#define EIKONALXX_GEOMETRY3D_HPP
#include <memory>
namespace EikonalXX
{
/// @class Geometry3D "geometry3d.hpp" "eikonalxx/geometry3d.hpp"
/// @brief Defines the 3D physical model geometry.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Geometry3D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Geometry3D();
    /// @brief Copy constructor.
    /// @param[in] geometry  The geometry class from which to initialize this
    ///                      class.
    Geometry3D(const Geometry3D &geometry);
    /// @brief Move constructor.
    /// @param[in,out] geometry  The geometry class from which to initialize
    ///                          this class.  On exit, geometry's behavior
    ///                          is undefined.
    Geometry3D(Geometry3D &&geometry) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] geometry   The geometry class to copy to this.
    /// @result A deep copy of the geometry class.
    Geometry3D& operator=(const Geometry3D &geometry);
    /// @brief Move assignment operator.
    /// @param[in,out] geometry   The geometry class whose memory will be moved
    ///                           to this.  On exit, geometry's behavior is
    ///                           undefined.
    /// @result The memory from geometry moved to this.
    Geometry3D& operator=(Geometry3D &&geometry) noexcept; 
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

    /// @name Properties in y direction
    /// @{
    /// @brief Sets the number of grid points in y.
    /// @param[in] ny  The number of grid points in y.  This must be at least 3.
    /// @throws std::invalid_argument if ny is too small.
    /// @note This is required by the solver.
    void setNumberOfGridPointsInY(int ny);
    /// @result The number of grid points in y.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInY() is false.
    [[nodiscard]] int getNumberOfGridPointsInY() const;
    /// @result The number of cells in x.
    /// @throw std::runtime_error if \c haveNumberOfGridPointsInY() is false.
    [[nodiscard]] int getNumberOfCellsInY() const;
    /// @result True indicates that the number of grid points in y is set.
    [[nodiscard]] bool haveNumberOfGridPointsInY() const noexcept;

    /// @brief Sets the grid spacing in y.  This is required for the solver.
    /// @param[in] dy   The grid spacing in y in meters.
    /// @throws std::invalid_argument if dy is not positive.
    /// @note This is required by the solver. 
    void setGridSpacingInY(double dy);
    /// @result The grid spacing in y in meters.
    /// @throws std::runtime_error if \c haveGridSpacingInY() is false.
    [[nodiscard]] double getGridSpacingInY() const;
    /// @result True indicates that the grid spacing in y was set. 
    [[nodiscard]] bool haveGridSpacingInY() const noexcept;

    /// @brief Sets the model's y origin.
    /// @param[in] y0   The model's y origin in meters.
    /// @note By default this is 0.
    void setOriginInY(double y0) noexcept;
    /// @result The model's y origin in meters.
    double getOriginInY() const noexcept; 
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
    /// @throws std::runtime_error if \c haveNumberOfGridPointsInX() is false,
    ///         \c haveNumberOfGridPointsInY() is false or,
    ///         \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfGridPoints() const;
    /// @result The number of cells in the geometry.
    /// @throws std::runtime_error if \c haveNumberOfGridPointsInX() is false,
    ///         \c haveNumberOfGridPointsInY() is false, or
    ///         \c haveNumberOfGridPointsInZ() is false.
    [[nodiscard]] int getNumberOfCells() const;

    /// @name Destructors
    /// @{
    /// @brief Destructor.
    ~Geometry3D();
    /// @brief Resets the class.
    void clear() noexcept;
    /// @} 
private:
    class Geometry3DImpl;
    std::unique_ptr<Geometry3DImpl> pImpl;
};
/// @result True if the left-hand side geometry is equivalent to the
///         right-hand side geometry.
bool operator==(const Geometry3D &lhs, const Geometry3D &rhs);
/// @result True if the left-hand side geometry is not equivalent to
///         the right-hand side geometry.
bool operator!=(const Geometry3D &lhs, const Geometry3D &rhs);
}
#endif
