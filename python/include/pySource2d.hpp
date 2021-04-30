#ifndef PYEIKONALXX_SOURCE2D_HPP
#define PYEIKONALXX_SOURCE2D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
class Source2D;
}
namespace PEikonalXX
{
class Geometry2D;
class Source2D
{
public:
    /// @brief Constructor.
    Source2D();
    /// @brief Copy c'tor.
    Source2D(const Source2D &source);
    Source2D(const EikonalXX::Source2D &source);
    /// @brief Copy c'tor.
    Source2D(Source2D &&source) noexcept;
    /// @brief Copy assignment.
    Source2D& operator=(const Source2D &source);
    Source2D& operator=(const EikonalXX::Source2D &source);
    /// @brief Move assignment.
    Source2D& operator=(Source2D &&source) noexcept;
    /// @brief Destructor.
    ~Source2D();
    /// @brief Reset the class.
    void clear() noexcept;

    void setGeometry(const Geometry2D &geometry);
    [[nodiscard]] bool haveGeometry() const noexcept;

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
    std::unique_ptr<EikonalXX::Source2D> pImpl;
};
void initializeSource2D(pybind11::module &m);
}
#endif
