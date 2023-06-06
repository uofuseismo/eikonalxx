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

/// Nodal velocities
void Model2D::setNodalVelocities(
    const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &velocities, 
    const EikonalXX::Ordering2D ordering)
{
    pybind11::buffer_info velocityBuffer = velocities.request();
    auto nGrid = static_cast<int> (velocityBuffer.size);
    const double *velocityPointer = (double *) (velocityBuffer.ptr);
    pImpl->setNodalVelocities(nGrid, velocityPointer, ordering); 
}

void Model2D::setCellularVelocities(
    const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> &velocities, 
    const EikonalXX::Ordering2D ordering)
{
    pybind11::buffer_info velocityBuffer = velocities.request();
    auto nCell = static_cast<int> (velocityBuffer.size);
    const double *velocityPointer = (double *) (velocityBuffer.ptr);
    pImpl->setCellularVelocities(nCell, velocityPointer, ordering); 
}

/// Write vtk
void Model2D::writeVTK(const std::string &fileName, const std::string &title)
{
    pImpl->writeVTK(fileName, title); 
}

/// Initialize class
void PEikonalXX::initializeModel2D(pybind11::module &module)
{
    pybind11::class_<PEikonalXX::Model2D> m(module, "VelocityModel2D");
    m.def(pybind11::init<> ());
    m.doc() = R""""(
This defines a 2D Cartesian velocity model.

Read-Only Properties
--------------------
is_initialized : bool
   True indicates the class is initialized.
geometry : Geometry2D
   The two-dimensionsal geometry.
)"""";

    m.def("__copy__", [](const Model2D &self)
    {
        return Model2D(self);
    });
    m.def("initialize",
          &Model2D::initialize, 
          "Initializes the class.");
    m.def_property_readonly("is_initialized",
                            &Model2D::isInitialized);
    m.def_property_readonly("geometry",
                            &Model2D::getGeometry);
    m.def("set_cellular_velocities",
          &Model2D::setCellularVelocities,
          "Sets the cell-based velocities.  The velocities are in m/s.",
          pybind11::arg("velocities"),
          pybind11::arg("ordering") = EikonalXX::Ordering2D::Natural);
    m.def("set_nodal_velocities",
          &Model2D::setNodalVelocities,
          "Sets the node-based velocities.  The velocities are in m/s.",
          pybind11::arg("velocities"),
          pybind11::arg("ordering") = EikonalXX::Ordering2D::Natural);
    m.def("write_vtk",
          &Model2D::writeVTK,
          "Writes the velocity model to a VTK file.",
          pybind11::arg("file_name"),
          pybind11::arg("title") = "velocity_m/s"); 
}
