#include "include/pyGeometry3d.hpp"
#include "eikonalxx/geometry3d.hpp"

using namespace PEikonalXX;

/// C'tor
Geometry3D::Geometry3D() :
    pImpl(std::make_unique<EikonalXX::Geometry3D> ())
{
}

/// Destructor
Geometry3D::~Geometry3D() = default;

/// Copy c'tor
Geometry3D::Geometry3D(const Geometry3D &geometry)
{
    *this = geometry;
}

Geometry3D::Geometry3D(const EikonalXX::Geometry3D &geometry)
{
    *this = geometry;
}

/// Move c'tor
Geometry3D::Geometry3D(Geometry3D &&geometry) noexcept
{
    *this = std::move(geometry);
}

/// Copy assignment 
Geometry3D& Geometry3D::operator=(const Geometry3D &geometry)
{
    if (&geometry == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Geometry3D> (*geometry.pImpl);
    return *this;
}

Geometry3D& Geometry3D::operator=(const EikonalXX::Geometry3D &geometry)
{
    pImpl = std::make_unique<EikonalXX::Geometry3D> (geometry);
    return *this;
}

/// Move assignment
Geometry3D& Geometry3D::operator=(Geometry3D &&geometry) noexcept
{
    if (&geometry == this){return *this;}
    pImpl = std::move(geometry.pImpl);
    return *this;
}

/// Clears the class
void Geometry3D::clear() noexcept
{
    pImpl->clear();
}

/// Grid spacing
void Geometry3D::setGridSpacingInX(const double dx)
{
    pImpl->setGridSpacingInX(dx);
}

double Geometry3D::getGridSpacingInX() const
{
    return pImpl->getGridSpacingInX();
}

void Geometry3D::setGridSpacingInY(const double dy)
{
    pImpl->setGridSpacingInY(dy);
}

double Geometry3D::getGridSpacingInY() const
{
    return pImpl->getGridSpacingInY();
}

void Geometry3D::setGridSpacingInZ(const double dz)
{
    pImpl->setGridSpacingInZ(dz);
}

double Geometry3D::getGridSpacingInZ() const
{
    return pImpl->getGridSpacingInZ();
}

/// Grid origin
void Geometry3D::setOriginInX(const double x0) noexcept
{
    pImpl->setOriginInX(x0);
}

double Geometry3D::getOriginInX() const noexcept
{
    return pImpl->getOriginInX();
}

void Geometry3D::setOriginInY(const double y0) noexcept
{
    pImpl->setOriginInY(y0);
}

double Geometry3D::getOriginInY() const noexcept
{
    return pImpl->getOriginInY();
}

void Geometry3D::setOriginInZ(const double z0) noexcept
{
    pImpl->setOriginInZ(z0);
}

double Geometry3D::getOriginInZ() const noexcept
{
    return pImpl->getOriginInZ();
}

/// Number of grid points
void Geometry3D::setNumberOfGridPointsInX(const int nx)
{
    pImpl->setNumberOfGridPointsInX(nx);
}

int Geometry3D::getNumberOfGridPointsInX() const
{
    return pImpl->getNumberOfGridPointsInX();
}

void Geometry3D::setNumberOfGridPointsInY(const int ny) 
{
    pImpl->setNumberOfGridPointsInY(ny);
}

int Geometry3D::getNumberOfGridPointsInY() const
{
    return pImpl->getNumberOfGridPointsInY();
}

void Geometry3D::setNumberOfGridPointsInZ(const int nz)
{
    pImpl->setNumberOfGridPointsInZ(nz);
}

int Geometry3D::getNumberOfGridPointsInZ() const
{
    return pImpl->getNumberOfGridPointsInZ();
}

/// Initialize class
void PEikonalXX::initializeGeometry3D(pybind11::module &m)
{
    pybind11::class_<PEikonalXX::Geometry3D> g(m, "Geometry3D");
    g.def(pybind11::init<> ());
    g.doc() = "This defines a 3D Cartesian grid-based geometry.  The geometry is in a left-handed coordinate system where increases +x right, +y increases away, and +z increases down.  There are three things that define a grid geometry - the number of grid points, the grid spacing, and, optionally, the grid origin.\n\nThe number of grid number of grid points is specified by (nx,ny,nz) where nx, ny, and nz must all be at least 3.\n\nThe grid spacing is given by (dx,dy,dz) in meters where dx, dy, and dz must be positive.\n\nThe origin is defined by (x0,y0,z0) in meters.  This will default to (0,0,0).";
    

    g.def_property("nx",
                   &Geometry3D::getNumberOfGridPointsInX,
                   &Geometry3D::setNumberOfGridPointsInX,
                   "The number of grid points in x.  This must be at least 3.");
    g.def_property("ny",
                   &Geometry3D::getNumberOfGridPointsInY,
                   &Geometry3D::setNumberOfGridPointsInY,
                   "The number of grid points in y.  This must be at least 3.");
    g.def_property("nz",
                   &Geometry3D::getNumberOfGridPointsInZ,
                   &Geometry3D::setNumberOfGridPointsInZ,
                   "The number of grid points in z.  This must be at least 3.");

    g.def_property("dx",
                   &Geometry3D::getGridSpacingInX,
                   &Geometry3D::setGridSpacingInX,
                   "The grid spacing in x in meters.  This must be positive.");
    g.def_property("dy",
                   &Geometry3D::getGridSpacingInY,
                   &Geometry3D::setGridSpacingInY,
                   "The grid spacing in y in meters.  This must be positive.");
    g.def_property("dz",
                   &Geometry3D::getGridSpacingInZ,
                   &Geometry3D::setGridSpacingInZ,
                   "The grid spacing in z in meters.  This must be positive.");

    g.def_property("x0",
                   &Geometry3D::getOriginInX,
                   &Geometry3D::setOriginInX,
                   "The grid origin in x in meters.");
    g.def_property("y0",
                   &Geometry3D::getOriginInY,
                   &Geometry3D::setOriginInY,
                   "The grid origin in y in meters.");
    g.def_property("z0",
                   &Geometry3D::getOriginInZ,
                   &Geometry3D::setOriginInZ,
                   "The grid origin in z in meters.");

    g.def("clear",
          &Geometry3D::clear,
          "Resets the class.");

    /// Pickling rules (makes this class copyable)
    g.def(pybind11::pickle(
        [](const Geometry3D &g) { //__getstate__
            auto x0 = g.getOriginInX();
            auto y0 = g.getOriginInY();
            auto z0 = g.getOriginInZ();
            auto dx = g.getGridSpacingInX(); // throws
            auto dy = g.getGridSpacingInY(); // throws
            auto dz = g.getGridSpacingInZ(); // throws
            auto nx = g.getNumberOfGridPointsInX(); // throws
            auto ny = g.getNumberOfGridPointsInY(); // throws
            auto nz = g.getNumberOfGridPointsInZ(); // throws
            return pybind11::make_tuple(x0, y0, z0, dx, dy, dz, nx, ny, nz);
        },
        [](pybind11::tuple t) { //__setstate__
            if (t.size() != 9)
            {
                throw std::runtime_error("Invalid state!");
            }
            auto x0 = t[0].cast<double> ();
            auto y0 = t[1].cast<double> ();
            auto z0 = t[2].cast<double> ();
            auto dx = t[3].cast<double> ();
            auto dy = t[4].cast<double> ();
            auto dz = t[5].cast<double> ();
            auto nx = t[6].cast<int> ();
            auto ny = t[7].cast<int> ();
            auto nz = t[8].cast<int> ();
            Geometry3D p;
            p.setOriginInX(x0);
            p.setOriginInY(y0);
            p.setOriginInZ(z0);
            p.setGridSpacingInX(dx);
            p.setGridSpacingInY(dy);
            p.setGridSpacingInZ(dz);
            p.setNumberOfGridPointsInX(nx);
            p.setNumberOfGridPointsInY(ny);
            p.setNumberOfGridPointsInZ(nz);
            return p;
        }
    ));
}
