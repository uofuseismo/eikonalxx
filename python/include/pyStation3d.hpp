#ifndef PYEIKONALXX_STATION_3D_HPP
#define PYEIKONALXX_STATION_3D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Station3D;
}
namespace PEikonalXX
{
class Geometry3D;
class Station3D
{
public:
    /// @brief Constructor.
    Station3D();
    /// @brief Copy c'tor.
    Station3D(const Station3D &station);
    Station3D(const EikonalXX::Station3D &station);
    /// @brief Copy c'tor.
    Station3D(Station3D &&station) noexcept;
    /// @brief Copy assignment.
    Station3D& operator=(const Station3D &station);
    Station3D& operator=(const EikonalXX::Station3D &station);
    /// @brief Move assignment.
    Station3D& operator=(Station3D &&station) noexcept;
    /// @brief Destructor.
    ~Station3D();
    /// @brief Reset the class.
    void clear() noexcept;

    void setGeometry(const Geometry3D &geometry);
    [[nodiscard]] bool haveGeometry() const noexcept;
    Geometry3D getGeometry() const;

    /// @brief X location
    void setLocationInX(double x);
    [[nodiscard]] bool haveLocationInX() const noexcept;
    [[nodiscard]] double getLocationInX() const;

    /// @brief Y location
    void setLocationInY(double y); 
    [[nodiscard]] bool haveLocationInY() const noexcept;
    [[nodiscard]] double getLocationInY() const;

    /// @brief Z location
    void setLocationInZ(double x);
    void setZToFreeSurface();
    [[nodiscard]] bool haveLocationInZ() const noexcept;
    [[nodiscard]] double getLocationInZ() const;

private:
    std::unique_ptr<EikonalXX::Station3D> pImpl;
};
void initializeStation3D(pybind11::module &m);
}
#endif
