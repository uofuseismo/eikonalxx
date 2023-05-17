#ifndef EIKONALXX_RAY_POINT_3D_HPP
#define EIKONALXX_RAY_POINT_3D_HPP
#include <memory>
namespace EikonalXX::Ray
{
/// @class Point3D "point3d.hpp" "eikonalxx/ray/point3d.hpp"
/// @brief Defines a point along a ray path.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Point3D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Point3D();
    /// @brief Copy constructor.
    /// @param[in] point  The point from which to initialize this class.
    Point3D(const Point3D &point);
    /// @brief Move constructor.
    /// @param[in,out] point  The point from which to initialize this class. 
    ///                       On exit, point's behavior is undefined.
    Point3D(Point3D &&point) noexcept;
    /// @brief Constructs a point from an (x, y, z) position.
    /// @param[in] x  The position in the model in x.  This has units of meters.
    /// @param[in] y  The position in the model in y.  This has units of meters.
    /// @param[in] z  The position in the model in z.  This has units of meters.
    Point3D(double x, double y, double z); 
    /// @}

    /// @name X Position
    /// @{

    /// @param[in] x  The position in the model in x.  This has units of meters.
    void setPositionInX(double x) noexcept;
    /// @result The position in the model in x (meters).
    [[nodiscard]] double getPositionInX() const;
    /// @result True indicates the x position was set.
    [[nodiscard]] bool havePositionInX() const noexcept;
    /// @}

    /// @name Y Position
    /// @{

    /// @param[in] y  The position in the model in y.  This has units of meters.
    void setPositionInY(double y) noexcept;
    /// @result The position in the model in y (meters).
    [[nodiscard]] double getPositionInY() const;
    /// @result True indicates the y position was set.
    [[nodiscard]] bool havePositionInY() const noexcept;
    /// @

    /// @name Z Position
    /// @{

    /// @param[in] z  The position in the model in z.  This has units of meters.
    void setPositionInZ(double z) noexcept;
    /// @result The position in the model in z (meters).
    [[nodiscard]] double getPositionInZ() const;
    /// @result True indicates the z position was set.
    [[nodiscard]] bool havePositionInZ() const noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] point  The point to copy to this.
    /// @result A deep copy of the input point.
    Point3D& operator=(const Point3D &point);
    /// @brief Move assignment operator.
    /// @param[in,out] point  The point whose memory will be moved to this.
    ///                       On exit, point's behavior is undefined.
    /// @result The memory from point moved to this.
    Point3D& operator=(Point3D &&point) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Point3D();
    /// @}
private:
    class Point3DImpl;
    std::unique_ptr<Point3DImpl> pImpl;
};
}
#endif 
