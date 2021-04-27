#ifndef PYEIKONALXX_GEOMETRY3D_HPP
#define PYEIKONALXX_GEOMETRY3D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Geometry3D;
}
namespace PEikonalXX
{
class Geometry3D
{
public:
    /// @brief Constructor.
    Geometry3D();
    /// @brief Copy c'tor.
    Geometry3D(const Geometry3D &geometry);
    Geometry3D(const EikonalXX::Geometry3D &geometry);
    /// @brief Copy c'tor.
    Geometry3D(Geometry3D &&geometry) noexcept;
    /// @brief Copy assignment.
    Geometry3D& operator=(const Geometry3D &geometry);
    Geometry3D& operator=(const EikonalXX::Geometry3D &geometry);
    /// @brief Move assignment.
    Geometry3D& operator=(Geometry3D &&geometry) noexcept;

    /// @brief Destructor
    ~Geometry3D();
    /// @brief Resets the class
    void clear() noexcept;

    /// @result Pointer to the native class.
    const EikonalXX::Geometry3D *getNativeClassPointer() const;
    /// Grid points in (x,y,z)
    void setNumberOfGridPointsInX(int nx);
    [[nodiscard]] int getNumberOfGridPointsInX() const;
    void setNumberOfGridPointsInY(int ny);
    [[nodiscard]] int getNumberOfGridPointsInY() const;
    void setNumberOfGridPointsInZ(int nz);
    [[nodiscard]] int getNumberOfGridPointsInZ() const;

     /// Grid spacing in (x,y,z)
     void setGridSpacingInX(double dx);
     [[nodiscard]] double getGridSpacingInX() const;
     void setGridSpacingInY(double dy);
     [[nodiscard]] double getGridSpacingInY() const;
     void setGridSpacingInZ(double dz);
     [[nodiscard]] double getGridSpacingInZ() const;

     // Grid origin in (x,y,z)
     void setOriginInX(double x0) noexcept;
     [[nodiscard]] double getOriginInX() const noexcept;
     void setOriginInY(double y0) noexcept;
     [[nodiscard]] double getOriginInY() const noexcept;
     void setOriginInZ(double z0) noexcept;
     [[nodiscard]] double getOriginInZ() const noexcept;
private:
    std::unique_ptr<EikonalXX::Geometry3D> pImpl;
};
void initializeGeometry3D(pybind11::module &m);
}
#endif
