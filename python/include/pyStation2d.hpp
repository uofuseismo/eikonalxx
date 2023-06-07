#ifndef PYEIKONALXX_STATION_2D_HPP
#define PYEIKONALXX_STATION_2D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Station2D;
}
namespace PEikonalXX
{
class Geometry2D;
class Station2D
{
public:
    /// @brief Constructor.
    Station2D();
    /// @brief Copy c'tor.
    Station2D(const Station2D &station);
    Station2D(const EikonalXX::Station2D &station);
    /// @brief Copy c'tor.
    Station2D(Station2D &&station) noexcept;
    /// @brief Copy assignment.
    Station2D& operator=(const Station2D &station);
    Station2D& operator=(const EikonalXX::Station2D &station);
    /// @brief Move assignment.
    Station2D& operator=(Station2D &&station) noexcept;
    /// @brief Destructor.
    ~Station2D();
    /// @brief Reset the class.
    void clear() noexcept;

    [[nodiscard]] const EikonalXX::Station2D *getNativeClassPointer() const;

    void setGeometry(const Geometry2D &geometry);
    [[nodiscard]] bool haveGeometry() const noexcept;
    Geometry2D getGeometry() const;

    /// @brief X location
    void setLocationInX(double x);
    [[nodiscard]] bool haveLocationInX() const noexcept;
    [[nodiscard]] double getLocationInX() const;

    /// @brief Z location
    void setLocationInZ(double x);
    void setZToFreeSurface();
    [[nodiscard]] bool haveLocationInZ() const noexcept;
    [[nodiscard]] double getLocationInZ() const;

private:
    std::unique_ptr<EikonalXX::Station2D> pImpl;
};
void initializeStation2D(pybind11::module &m);
}
#endif
