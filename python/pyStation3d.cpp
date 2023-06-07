#include "include/pyStation3d.hpp"
#include "include/pyGeometry3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/station3d.hpp"

using namespace PEikonalXX;

/// C'tor
Station3D::Station3D() :
    pImpl(std::make_unique<EikonalXX::Station3D> ())
{
}

Station3D::Station3D(const Station3D &station)
{
    *this = station;
}

Station3D::Station3D(const EikonalXX::Station3D &station)
{
    *this = station;
}

Station3D::Station3D(Station3D &&station) noexcept
{
    *this = std::move(station);
}

/// Destructor
Station3D::~Station3D() = default;

/// Operators
Station3D& Station3D::operator=(const Station3D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Station3D> (*source.pImpl);
    return *this;
}

Station3D& Station3D::operator=(const EikonalXX::Station3D &source)
{
    pImpl = std::make_unique<EikonalXX::Station3D> (source);
    return *this;
}

Station3D& Station3D::operator=(Station3D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

/// Set geometry
void Station3D::setGeometry(const Geometry3D &geometry)
{
    auto geoNative = geometry.getNativeClassPointer();
    pImpl->setGeometry(*geoNative);
}

bool Station3D::haveGeometry() const noexcept
{
    return pImpl->haveGeometry();
}

Geometry3D Station3D::getGeometry() const
{
    if (!haveGeometry())
    {
        throw std::invalid_argument("Geometry not yet set");
    }
    Geometry3D geometry(pImpl->getGeometry());
    return geometry;
}

/// X location
void Station3D::setLocationInX(const double x)
{
    pImpl->setLocationInX(x);
}

bool Station3D::haveLocationInX() const noexcept
{
    return pImpl->haveLocationInX();
}

double Station3D::getLocationInX() const
{
    return pImpl->getLocationInX();
}

/// Y location
void Station3D::setLocationInY(const double y)
{
    pImpl->setLocationInY(y);
}

bool Station3D::haveLocationInY() const noexcept
{
    return pImpl->haveLocationInY();
}

double Station3D::getLocationInY() const
{
    return pImpl->getLocationInY();
}

/// Z location
void Station3D::setLocationInZ(const double z)
{
    pImpl->setLocationInZ(z);
}

void Station3D::setZToFreeSurface()
{
    pImpl->setZToFreeSurface();
}

bool Station3D::haveLocationInZ() const noexcept
{
    return pImpl->haveLocationInZ();
}

double Station3D::getLocationInZ() const
{
    return pImpl->getLocationInZ();
}

const EikonalXX::Station3D* Station3D::getNativeClassPointer() const
{
    return pImpl.get();
}

/// Reset
void Station3D::clear() noexcept
{
    pImpl->clear();
}

/// Initialize class
void PEikonalXX::initializeStation3D(pybind11::module &m) 
{
    pybind11::class_<PEikonalXX::Station3D> s(m, "Station3D");
    s.def(pybind11::init<> ());
    s.doc() = R"""(
This defines a station in a 3D Cartesian grid-based geometry.
To use this you must first set the geometry.

Properties
----------
geometry : Geometry3D
    The 3D model geometry.
x : float
    The station's position (meters) in x in the geometry.
y : float
    The station's position (meters) in y in the geometry.
z : float
    The station's position (meters) in z in the geometry.
    Note, this can be forced to the free-surface of the geometry
    with set_z_to_free_surface().
)""";
    s.def("__copy__", [](const Station3D &self)
    {
        return Station3D(self);
    });
    s.def_property("geometry",
                   &Station3D::getGeometry,
                   &Station3D::setGeometry);
    s.def_property("x",
                   &Station3D::getLocationInX,
                   &Station3D::setLocationInX);
    s.def_property("y",
                   &Station3D::getLocationInY,
                   &Station3D::setLocationInY);
    s.def_property("z",
                   &Station3D::getLocationInZ,
                   &Station3D::setLocationInZ);

    s.def("set_z_to_free_surface",
          &Station3D::setZToFreeSurface,
          "Sets the station's position in z to the free surface.");
    s.def("clear",
          &Station3D::clear,
          "Resets the class.");
    /// Pickling rules
    s.def(pybind11::pickle(
        [](const Station3D &g) { //__getstate__
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
            Station3D s;
            s.setGeometry(geometry);
            s.setLocationInX(x);
            s.setLocationInZ(z);
            return s;
        }
    ));
}
