#include "include/pyGeometry2d.hpp"
#include "eikonalxx/geometry2d.hpp"

using namespace PEikonalXX;

/// C'tor
Geometry2D::Geometry2D() :
    pImpl(std::make_unique<EikonalXX::Geometry2D> ())
{
}

/// Destructor
Geometry2D::~Geometry2D() = default;

/// Copy c'tor
Geometry2D::Geometry2D(const Geometry2D &geometry)
{
    *this = geometry;
}

Geometry2D::Geometry2D(const EikonalXX::Geometry2D &geometry)
{
    *this = geometry;
}

/// Move c'tor
Geometry2D::Geometry2D(Geometry2D &&geometry) noexcept
{
    *this = std::move(geometry);
}

/// Copy assignment 
Geometry2D& Geometry2D::operator=(const Geometry2D &geometry)
{
    if (&geometry == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Geometry2D> (*geometry.pImpl);
    return *this;
}

Geometry2D& Geometry2D::operator=(const EikonalXX::Geometry2D &geometry)
{
    pImpl = std::make_unique<EikonalXX::Geometry2D> (geometry);
    return *this;
}

/// Move assignment
Geometry2D& Geometry2D::operator=(Geometry2D &&geometry) noexcept
{
    if (&geometry == this){return *this;}
    pImpl = std::move(geometry.pImpl);
    return *this;
}

/// Clears the class
void Geometry2D::clear() noexcept
{
    pImpl->clear();
}

/// Get pointer to the native class
const EikonalXX::Geometry2D* Geometry2D::getNativeClassPointer() const
{
    return pImpl.get();
}

/// Grid spacing
void Geometry2D::setGridSpacingInX(const double dx)
{
    pImpl->setGridSpacingInX(dx);
}

double Geometry2D::getGridSpacingInX() const
{
    return pImpl->getGridSpacingInX();
}

void Geometry2D::setGridSpacingInZ(const double dz)
{
    pImpl->setGridSpacingInZ(dz);
}

double Geometry2D::getGridSpacingInZ() const
{
    return pImpl->getGridSpacingInZ();
}

/// Grid origin
void Geometry2D::setOriginInX(const double x0) noexcept
{
    pImpl->setOriginInX(x0);
}

double Geometry2D::getOriginInX() const noexcept
{
    return pImpl->getOriginInX();
}

void Geometry2D::setOriginInZ(const double z0) noexcept
{
    pImpl->setOriginInZ(z0);
}

double Geometry2D::getOriginInZ() const noexcept
{
    return pImpl->getOriginInZ();
}

/// Number of grid points
void Geometry2D::setNumberOfGridPointsInX(const int nx)
{
    pImpl->setNumberOfGridPointsInX(nx);
}

int Geometry2D::getNumberOfGridPointsInX() const
{
    return pImpl->getNumberOfGridPointsInX();
}

void Geometry2D::setNumberOfGridPointsInZ(const int nz)
{
    pImpl->setNumberOfGridPointsInZ(nz);
}

int Geometry2D::getNumberOfGridPointsInZ() const
{
    return pImpl->getNumberOfGridPointsInZ();
}

/// Initialize class
void PEikonalXX::initializeGeometry2D(pybind11::module &m)
{
    pybind11::class_<PEikonalXX::Geometry2D> g(m, "Geometry2D");
    g.def(pybind11::init<> ());
    g.doc() = R"""(
This defines a 2D Cartesian grid-based geometry.  The geometry is increases +x right and +z down. 
There are three things that define a grid geometry - the number of grid points, the grid spacing, and, optionally, the grid origin.

Required Properties
-------------------
nx : int
    The number of grid points in x.  This must be at least 3.
nz : int
    The number of grid points in z.  This must be at least 3.
dx : float
    The grid spacing in meters in x.  This must be positive.
dz : float
    The grid spacing in meters in z.  This must be positive.

Optional Properties
-------------------
x0 : float
    The x origin in meters.  By default this is 0.
z0 : float
    The z origin in meters.  By default this is 0.

)""";
    g.def("__copy__", [](const Geometry2D &self)
    {   
        return Geometry2D(self);
    });
    g.def_property("nx",
                   &Geometry2D::getNumberOfGridPointsInX,
                   &Geometry2D::setNumberOfGridPointsInX);
    g.def_property("nz",
                   &Geometry2D::getNumberOfGridPointsInZ,
                   &Geometry2D::setNumberOfGridPointsInZ);

    g.def_property("dx",
                   &Geometry2D::getGridSpacingInX,
                   &Geometry2D::setGridSpacingInX);
    g.def_property("dz",
                   &Geometry2D::getGridSpacingInZ,
                   &Geometry2D::setGridSpacingInZ);

    g.def_property("x0",
                   &Geometry2D::getOriginInX,
                   &Geometry2D::setOriginInX);
    g.def_property("z0",
                   &Geometry2D::getOriginInZ,
                   &Geometry2D::setOriginInZ);

    g.def("clear",
          &Geometry2D::clear,
          "Resets the class.");

    /// Pickling rules (makes this class copyable)
    g.def(pybind11::pickle(
        [](const Geometry2D &g) { //__getstate__
            auto x0 = g.getOriginInX();
            auto z0 = g.getOriginInZ();
            auto dx = g.getGridSpacingInX(); // throws
            auto dz = g.getGridSpacingInZ(); // throws
            auto nx = g.getNumberOfGridPointsInX(); // throws
            auto nz = g.getNumberOfGridPointsInZ(); // throws
            return pybind11::make_tuple(x0, z0, dx, dz, nx, nz);
        },
        [](pybind11::tuple t) { //__setstate__
            if (t.size() != 6)
            {
                throw std::runtime_error("Invalid state!");
            }
            auto x0 = t[0].cast<double> ();
            auto z0 = t[1].cast<double> ();
            auto dx = t[2].cast<double> ();
            auto dz = t[3].cast<double> ();
            auto nx = t[4].cast<int> ();
            auto nz = t[5].cast<int> ();
            Geometry2D p;
            p.setOriginInX(x0);
            p.setOriginInZ(z0);
            p.setGridSpacingInX(dx);
            p.setGridSpacingInZ(dz);
            p.setNumberOfGridPointsInX(nx);
            p.setNumberOfGridPointsInZ(nz);
            return p;
        }
    ));
}
