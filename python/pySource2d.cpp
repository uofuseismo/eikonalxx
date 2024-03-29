#include "include/pySource2d.hpp"
#include "include/pyGeometry2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"

using namespace PEikonalXX;

/// C'tor
Source2D::Source2D() :
    pImpl(std::make_unique<EikonalXX::Source2D> ())
{
}

Source2D::Source2D(const Source2D &source)
{
    *this = source;
}

Source2D::Source2D(const EikonalXX::Source2D &source)
{
    *this = source;
}

Source2D::Source2D(Source2D &&source) noexcept
{
    *this = std::move(source);
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

const EikonalXX::Source2D* Source2D::getNativeClassPointer() const
{
    return pImpl.get();
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

Geometry2D Source2D::getGeometry() const
{
    if (!haveGeometry())
    {
        throw std::invalid_argument("Geometry not yet set");
    }
    Geometry2D geometry(pImpl->getGeometry());
    return geometry;
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
    s.doc() = R"""(
This defines a source in a 2D Cartesian grid-based geometry.  To use this you must first set the geometry.

Properties
----------
geometry : Geometry2D
    The 2D model geometry.
x : float
    The source's position (meters) in x in the geometry.
z : float
    The source's position (meters) in z in the geometry.
    Note, this can be forced to the free-surface of the geometry
    with source2d.set_z_to_free_surface().
)""";
    s.def("__copy__", [](const Source2D &self)
    {
        return Source2D(self);
    });
    s.def_property("geometry",
                   &Source2D::getGeometry,
                   &Source2D::setGeometry);
    s.def_property("x",
                   &Source2D::getLocationInX,
                   &Source2D::setLocationInX);
    s.def_property("z",
                   &Source2D::getLocationInZ,
                   &Source2D::setLocationInZ);
    s.def("set_z_to_free_surface",
          &Source2D::setZToFreeSurface,
          "Sets the source's position in z to the free surface.");
    s.def("clear",
          &Source2D::clear,
          "Resets the class.");

    /// Pickling rules
    s.def(pybind11::pickle(
        [](const Source2D &g) { //__getstate__
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
            auto geometry = t[2].cast<Geometry2D> ();
            Source2D s;
            s.setGeometry(geometry);
            s.setLocationInX(x);
            s.setLocationInZ(z);
            return s;
        }
    ));
}
