#include "include/pySource2d.hpp"
#include "include/pyGeometry2d.hpp"
#include "eikonalxx/source2d.hpp"

using namespace PEikonalXX;

/// C'tor
Source2D::Source2D() :
    pImpl(std::make_unique<EikonalXX::Source2D> ())
{
}

/// Destructor
Source2D::~Source2D() = default;

/// Operators
Source2D& Source2D::operator=(const Source2D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Source2D> (*source.pImpl);
    return *this;
}

Source2D& Source2D::operator=(const EikonalXX::Source2D &source)
{
    pImpl = std::make_unique<EikonalXX::Source2D> (source);
    return *this;
}

Source2D& Source2D::operator=(Source2D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

/// Set geometry
void Source2D::setGeometry(const Geometry2D &geometry)
{
    auto geoNative = geometry.getNativeClassPointer();
    pImpl->setGeometry(*geoNative);
}

bool Source2D::haveGeometry() const noexcept
{
    return pImpl->haveGeometry();
}

/// X location
void Source2D::setLocationInX(const double x)
{
    pImpl->setLocationInX(x);
}

bool Source2D::haveLocationInX() const noexcept
{
    return pImpl->haveLocationInX();
}

double Source2D::getLocationInX() const
{
    return pImpl->getLocationInX();
}

/// Z location
void Source2D::setLocationInZ(const double z)
{
    pImpl->setLocationInZ(z);
}

void Source2D::setZToFreeSurface()
{
    pImpl->setZToFreeSurface();
}

bool Source2D::haveLocationInZ() const noexcept
{
    return pImpl->haveLocationInZ();
}

double Source2D::getLocationInZ() const
{
    return pImpl->getLocationInZ();
}

/// Reset
void Source2D::clear() noexcept
{
    pImpl->clear();
}

/// Initialize class
void PEikonalXX::initializeSource2D(pybind11::module &m) 
{
    pybind11::class_<PEikonalXX::Source2D> s(m, "Source2D");
    s.def(pybind11::init<> ());
    s.doc() = "This defines a source in a 2D Cartesian grid-based geometry.  To use this you must first set the geometry.  Then the source's x position in the geometry (meters) can be set.  Likewise, the source's z position in the geometry (meters) can be set.  Furthermore, the source's z location can be fixed to the top of the model, i.e., the `free surface', using set_z_to_free_surface.";

    s.def("set_geometry",
          &Source2D::setGeometry,
          "Sets the model geometry.");

    s.def_property("x",
                   &Source2D::getLocationInX,
                   &Source2D::setLocationInX,
                   "The source's x position");
    s.def_property("z",
                   &Source2D::getLocationInZ,
                   &Source2D::setLocationInZ,
                   "The source's z position");
    s.def("set_z_to_free_surface",
          &Source2D::setZToFreeSurface,
          "Sets the source's position in z to the free surface.");

    s.def("clear",
          &Source2D::clear,
          "Resets the class.");
}
