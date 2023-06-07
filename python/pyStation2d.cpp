#include "include/pyStation2d.hpp"
#include "include/pyGeometry2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/station2d.hpp"

using namespace PEikonalXX;

/// C'tor
Station2D::Station2D() :
    pImpl(std::make_unique<EikonalXX::Station2D> ())
{
}

Station2D::Station2D(const Station2D &station)
{
    *this = station;
}

Station2D::Station2D(const EikonalXX::Station2D &station)
{
    *this = station;
}

Station2D::Station2D(Station2D &&station) noexcept
{
    *this = std::move(station);
}

/// Destructor
Station2D::~Station2D() = default;

/// Operators
Station2D& Station2D::operator=(const Station2D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Station2D> (*source.pImpl);
    return *this;
}

Station2D& Station2D::operator=(const EikonalXX::Station2D &source)
{
    pImpl = std::make_unique<EikonalXX::Station2D> (source);
    return *this;
}

Station2D& Station2D::operator=(Station2D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

const EikonalXX::Station2D* Station2D::getNativeClassPointer() const
{
    return pImpl.get();
}

/// Set geometry
void Station2D::setGeometry(const Geometry2D &geometry)
{
    auto geoNative = geometry.getNativeClassPointer();
    pImpl->setGeometry(*geoNative);
}

bool Station2D::haveGeometry() const noexcept
{
    return pImpl->haveGeometry();
}

Geometry2D Station2D::getGeometry() const
{
    if (!haveGeometry())
    {
        throw std::invalid_argument("Geometry not yet set");
    }
    Geometry2D geometry(pImpl->getGeometry());
    return geometry;
}

/// X location
void Station2D::setLocationInX(const double x)
{
    pImpl->setLocationInX(x);
}

bool Station2D::haveLocationInX() const noexcept
{
    return pImpl->haveLocationInX();
}

double Station2D::getLocationInX() const
{
    return pImpl->getLocationInX();
}

/// Z location
void Station2D::setLocationInZ(const double z)
{
    pImpl->setLocationInZ(z);
}

void Station2D::setZToFreeSurface()
{
    pImpl->setZToFreeSurface();
}

bool Station2D::haveLocationInZ() const noexcept
{
    return pImpl->haveLocationInZ();
}

double Station2D::getLocationInZ() const
{
    return pImpl->getLocationInZ();
}

/// Reset
void Station2D::clear() noexcept
{
    pImpl->clear();
}

/// Initialize class
void PEikonalXX::initializeStation2D(pybind11::module &m) 
{
    pybind11::class_<PEikonalXX::Station2D> s(m, "Station2D");
    s.def(pybind11::init<> ());
    s.doc() = R"""(
This defines a station in a 2D Cartesian grid-based geometry.
To use this you must first set the geometry.

Properties
----------
geometry : Geometry2D
    The 2D model geometry.
x : float
    The station's position (meters) in x in the geometry.
z : float
    The station's position (meters) in z in the geometry.
    Note, this can be forced to the free-surface of the geometry
    with set_z_to_free_surface().
)""";
    s.def("__copy__", [](const Station2D &self)
    {
        return Station2D(self);
    });
    s.def_property("geometry",
                   &Station2D::getGeometry,
                   &Station2D::setGeometry);
    s.def_property("x",
                   &Station2D::getLocationInX,
                   &Station2D::setLocationInX);
    s.def_property("z",
                   &Station2D::getLocationInZ,
                   &Station2D::setLocationInZ);

    s.def("set_z_to_free_surface",
          &Station2D::setZToFreeSurface,
          "Sets the station's position in z to the free surface.");
    s.def("clear",
          &Station2D::clear,
          "Resets the class.");
    /// Pickling rules
    s.def(pybind11::pickle(
        [](const Station2D &g) { //__getstate__
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
            Station2D s;
            s.setGeometry(geometry);
            s.setLocationInX(x);
            s.setLocationInZ(z);
            return s;
        }
    ));
}
