#ifndef PYEIKONALXX_SOURCE3D_HPP
#define PYEIKONALXX_SOURCE3D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Source3D;
}
namespace PEikonalXX
{
class Geometry3D;
class Source3D
{
public:
    /// @brief Constructor.
    Source3D();
    /// @brief Copy c'tor.
    Source3D(const Source3D &source);
    Source3D(const EikonalXX::Source3D &source);
    /// @brief Copy c'tor.
    Source3D(Source3D &&source) noexcept;
    /// @brief Copy assignment.
    Source3D& operator=(const Source3D &source);
    Source3D& operator=(const EikonalXX::Source3D &source);
    /// @brief Move assignment.
    Source3D& operator=(Source3D &&source) noexcept;
    /// @result Native class
    const EikonalXX::Source3D* getNativeClassPointer() const;
    /// @brief Destructor.
    ~Source3D();
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
    std::unique_ptr<EikonalXX::Source3D> pImpl;
};
void initializeSource3D(pybind11::module &m);
}
#endif
