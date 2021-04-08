#ifndef PYEIKONALXX_GEOMETRY2D_HPP
#define PYEIKONALXX_GEOMETRY2D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Geometry2D;
}
namespace PEikonalXX
{
class Geometry2D
{
public:
    /// @brief Constructor.
    Geometry2D();
    /// @brief Copy c'tor.
    Geometry2D(const Geometry2D &geometry);
    Geometry2D(const EikonalXX::Geometry2D &geometry);
    /// @brief Copy c'tor.
    Geometry2D(Geometry2D &&geometry) noexcept;
    /// @brief Copy assignment.
    Geometry2D& operator=(const Geometry2D &geometry);
    Geometry2D& operator=(const EikonalXX::Geometry2D &geometry);
    /// @brief Move assignment.
    Geometry2D& operator=(Geometry2D &&geometry) noexcept;

    /// @brief Destructor
    ~Geometry2D();
    /// @brief Resets the class
    void clear() noexcept;

    /// Grid points in (x,z)
    void setNumberOfGridPointsInX(int nx);
    [[nodiscard]] int getNumberOfGridPointsInX() const;
    void setNumberOfGridPointsInZ(int nz);
    [[nodiscard]] int getNumberOfGridPointsInZ() const;

     /// Grid spacing in (x,z)
     void setGridSpacingInX(double dx);
     [[nodiscard]] double getGridSpacingInX() const;
     void setGridSpacingInZ(double dz);
     [[nodiscard]] double getGridSpacingInZ() const;

     // Grid origin in (x,z)
     void setOriginInX(double x0) noexcept;
     [[nodiscard]] double getOriginInX() const noexcept;
     void setOriginInZ(double z0) noexcept;
     [[nodiscard]] double getOriginInZ() const noexcept;
private:
    std::unique_ptr<EikonalXX::Geometry2D> pImpl;
};
void initializeGeometry2D(pybind11::module &m);
}
#endif
