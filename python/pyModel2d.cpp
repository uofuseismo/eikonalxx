#include <string>
#include <vector>
#include <pybind11/stl.h>
#include "include/pyModel2d.hpp"
#include "include/pyGeometry2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/geometry2d.hpp"

using namespace PEikonalXX;

/// C'tor
Model2D::Model2D() :
    pImpl(std::make_unique<EikonalXX::Model2D<double>> ()) 
{
}

/// Destructor
Model2D::~Model2D() = default;

/// Reset class/release memory
void Model2D::clear() noexcept
{
    pImpl->clear();
}

/// Copy c'tor
Model2D::Model2D(const Model2D &model)
{
    *this = model;
}

/// Copy c'tor
Model2D::Model2D(const EikonalXX::Model2D<double> &model)
{
    *this = model;
}

/// Move c'tor
Model2D::Model2D(Model2D &&model) noexcept
{
    *this = std::move(model);
}

/// Copy
Model2D& Model2D::operator=(const Model2D &model)
{
    if (&model == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Model2D<double>> (*model.pImpl);
    return *this;
}

/// Copy
Model2D& Model2D::operator=(const EikonalXX::Model2D<double> &model)
{
    pImpl = std::make_unique<EikonalXX::Model2D<double>> (model);
    return *this;
}

/// Move assignment
Model2D& Model2D::operator=(Model2D &&model) noexcept
{
    if (&model == this){return *this;}
    pImpl = std::move(model.pImpl);
    return *this;
}

/// Get pointer to the native class
const EikonalXX::Model2D<double>* Model2D::getNativeClassPointer() const
{
    return pImpl.get();
}

/// Initialize
void Model2D::initialize(const Geometry2D &geometryIn)
{
    EikonalXX::Geometry2D geometry(*geometryIn.getNativeClassPointer());
    pImpl->initialize(geometry); 
}

/// Initialized?
bool Model2D::isInitialized() const noexcept
{
    return pImpl->isInitialized();
}

/// Get the model geometry
Geometry2D Model2D::getGeometry() const
{
    Geometry2D geometry(pImpl->getGeometry());
    return geometry;
}

/// Initialize class
void PEikonalXX::initializeModel2D(pybind11::module &module)
{
    pybind11::class_<PEikonalXX::Model2D> m(module, "VelocityModel2D");
    m.def(pybind11::init<> ());
    m.doc() = "This defines a 2D Cartesian velocity model.";

    m.def("initialize",
          &Model2D::initialize, 
          "Initializes the class.");
    m.def("is_initialized",
          &Model2D::isInitialized,
          "True indicates that the class is initialized");
    m.def("get_geometry",
          &Model2D::getGeometry,
          "Gets the underlying geometry.");

}
