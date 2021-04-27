#ifndef PYEIKONALXX_MODEL2D_HPP
#define PYEIKONALXX_MODEL2D_HPP
#include <memory>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
template<class T> class Model2D;
}
namespace PEikonalXX
{
class Geometry2D;
class Model2D
{
public:
    /// @brief Constructor.
    Model2D();
    /// @brief Copy c'tor.
    Model2D(const Model2D &model);
    Model2D(const EikonalXX::Model2D<double> &model);
    /// @brief Copy c'tor.
    Model2D(Model2D &&model) noexcept;
    /// @brief Copy assignment.
    Model2D& operator=(const Model2D &model);
    Model2D& operator=(const EikonalXX::Model2D<double> &model);
    /// @brief Move assignment.
    Model2D& operator=(Model2D &&model) noexcept;

    /// @result Pointer to the native class.
    const EikonalXX::Model2D<double> *getNativeClassPointer() const;
    /// @brief Destructor.
    ~Model2D();
    /// @brief Resets the class and releases memory.
    void clear() noexcept;

    /// @brief Initialize.
    void initialize(const Geometry2D &geometry);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @result The number of grid points.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfGridPoints() const; 
    /// @result The number of grid points.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfCells() const;
    /// @result The model geometry.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] Geometry2D getGeometry() const;
private:
    std::unique_ptr<EikonalXX::Model2D<double>> pImpl;
};
void initializeModel2D(pybind11::module &m);
}
#endif
