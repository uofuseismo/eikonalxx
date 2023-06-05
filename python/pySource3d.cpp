#include "include/pySource3d.hpp"
#include "include/pyGeometry3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/source3d.hpp"

using namespace PEikonalXX;

/// C'tor
Source3D::Source3D() :
    pImpl(std::make_unique<EikonalXX::Source3D> ())
{
}

Source3D::Source3D(const Source3D &source)
{
    *this = source;
}

Source3D::Source3D(const EikonalXX::Source3D &source)
{
    *this = source;
}

Source3D::Source3D(Source3D &&source) noexcept
{
    *this = std::move(source);
}

/// Destructor
Source3D::~Source3D() = default;

/// Operators
Source3D& Source3D::operator=(const Source3D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Source3D> (*source.pImpl);
    return *this;
}

Source3D& Source3D::operator=(const EikonalXX::Source3D &source)
{
    pImpl = std::make_unique<EikonalXX::Source3D> (source);
    return *this;
}

Source3D& Source3D::operator=(Source3D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

const EikonalXX::Source3D* Source3D::getNativeClassPointer() const
{
    return pImpl.get();
}

/// Set geometry
void Source3D::setGeometry(const Geometry3D &geometry)
{
    auto geoNative = geometry.getNativeClassPointer();
    pImpl->setGeometry(*geoNative);
}

bool Source3D::haveGeometry() const noexcept
{
    return pImpl->haveGeometry();
}

Geometry3D Source3D::getGeometry() const
{
    if (!haveGeometry())
    {
        throw std::invalid_argument("Geometry not yet set");
    }
    Geometry3D geometry(pImpl->getGeometry());
    return geometry;
}

/// X location
void Source3D::setLocationInX(const double x)
{
    pImpl->setLocationInX(x);
}

bool Source3D::haveLocationInX() const noexcept
{
    return pImpl->haveLocationInX();
}

double Source3D::getLocationInX() const
{
    return pImpl->getLocationInX();
}

/// Y location
void Source3D::setLocationInY(const double y)
{
    pImpl->setLocationInY(y);
}

bool Source3D::haveLocationInY() const noexcept
{
    return pImpl->haveLocationInY();
}

double Source3D::getLocationInY() const
{
    return pImpl->getLocationInY();
}

/// Z location
void Source3D::setLocationInZ(const double z)
{
    pImpl->setLocationInZ(z);
}

void Source3D::setZToFreeSurface()
{
    pImpl->setZToFreeSurface();
}

bool Source3D::haveLocationInZ() const noexcept
{
    return pImpl->haveLocationInZ();
}

double Source3D::getLocationInZ() const
{
    return pImpl->getLocationInZ();
}

/// Reset
void Source3D::clear() noexcept
{
    pImpl->clear();
}

/// Initialize class
void PEikonalXX::initializeSource3D(pybind11::module &m) 
{
    pybind11::class_<PEikonalXX::Source3D> s(m, "Source3D");
    s.def(pybind11::init<> ());
    s.doc() = R"""(
This defines a source in a 3D Cartesian grid-based geometry.  To use this you must first set the geometry.

Properties
----------
geometry : Geometry3D
    The 3D model geometry.
x : float
    The source's position (meters) in x in the geometry.
y : float
    The source's position (meters) in y in the geometry.
z : float
    The source's position (meters) in z in the geometry.
    Note, this can be forced to the free-surface of the geometry
    with source2d.set_z_to_free_surface().
)""";
    s.def("__copy__", [](const Source3D &self)
    {
        return Source3D(self);
    });
    s.def_property("geometry",
                   &Source3D::getGeometry,
                   &Source3D::setGeometry);
    s.def_property("x",
                   &Source3D::getLocationInX,
                   &Source3D::setLocationInX);
    s.def_property("y",
                   &Source3D::getLocationInY,
                   &Source3D::setLocationInY);
    s.def_property("z",
                   &Source3D::getLocationInZ,
                   &Source3D::setLocationInZ);
    s.def("set_z_to_free_surface",
          &Source3D::setZToFreeSurface,
          "Sets the source's position in z to the free surface.");
    s.def("clear",
          &Source3D::clear,
          "Resets the class.");

    /// Pickling rules
    s.def(pybind11::pickle(
        [](const Source3D &g) { //__getstate__
            auto x = g.getLocationInX();
            auto z = g.getLocationInZ();
            auto geometry = g.getGeometry();
            return pybind11::make_tuple(x, z, geometry);
        },
        [](pybind11::tuple t) { //__setstate__
            if (t.size() != 3)
            {
                throw std::runtime_error("Invalid state!");
            }
            auto x = t[0].cast<double> ();
            auto z = t[1].cast<double> ();
            auto geometry = t[2].cast<Geometry3D> ();
            Source3D s;
            s.setGeometry(geometry);
            s.setLocationInX(x);
            s.setLocationInZ(z);
            return s;
        }
    ));
}
