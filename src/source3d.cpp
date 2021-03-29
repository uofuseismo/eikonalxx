#include "eikonalxx/source3d.hpp"
#include "eikonalxx/geometry3d.hpp"

using namespace EikonalXX;

class Source3D::Source3DImpl
{
public:
    Geometry3D mGeometry;
    double mX = 0; 
    double mY = 0;
    double mZ = 0;
    double mXOffset = 0;
    double mYOffset = 0;
    double mZOffset = 0;
    int mCellX = 0;
    int mCellY = 0;
    int mCellZ = 0;
    bool mHaveXLocation = false;
    bool mHaveYLocation = false;
    bool mHaveZLocation = false;
    bool mHaveGeometry = false;
};

/// Reset the class
void Source3D::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mX = 0;
    pImpl->mY = 0;
    pImpl->mZ = 0;
    pImpl->mXOffset = 0;
    pImpl->mYOffset = 0;
    pImpl->mZOffset = 0;
    pImpl->mCellX = 0;
    pImpl->mCellY = 0;
    pImpl->mCellZ = 0;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveYLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = false;
}

/// C'tor
Source3D::Source3D() :
    pImpl(std::make_unique<Source3DImpl> ())
{
}

/// Copy c'tor
Source3D::Source3D(const Source3D &source)
{
    *this = source;
}

/// Move c'tor
Source3D::Source3D(Source3D &&source) noexcept
{
    *this = std::move(source);
}

/// Copy assignment operator
Source3D& Source3D::operator=(const Source3D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<Source3DImpl> (*source.pImpl);
    return *this;
}

/// Move assignment operator
Source3D& Source3D::operator=(Source3D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

/// Destructor
Source3D::~Source3D() = default;

/// Geometry
void Source3D::setGeometry(const Geometry3D &geometry)
{
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x must be set");
    }
    if (!geometry.haveGridSpacingInY())
    {
        throw std::invalid_argument("Grid spacing in y must be set");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("Grid spacing in z must be set");
    }
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Number of grid points in x must be set");
    }
    if (!geometry.haveNumberOfGridPointsInY())
    {
        throw std::invalid_argument("Number of grid points in y must be set");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Number of grid points in z must be set");
    }
    pImpl->mGeometry = geometry;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveYLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = true;
}

Geometry3D Source3D::getGeometry() const
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    return pImpl->mGeometry;
}

bool Source3D::haveGeometry() const noexcept
{
    return pImpl->mHaveGeometry;
}

/// x location
void Source3D::setLocationInX(const double x)
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    double dx = pImpl->mGeometry.getGridSpacingInX();
    double x0 = pImpl->mGeometry.getOriginInX();
    double x1 = x0 + static_cast<double>  (nx - 1)*dx;
    if (x < x0 || x > x1)
    {
        throw std::invalid_argument("x = " + std::to_string(x)
                                  + " must be in range [" + std::to_string(x0)
                                  + "," + std::to_string(x1) + "]");
    }
    pImpl->mX = x;
    pImpl->mXOffset = x - x0;
    pImpl->mCellX = static_cast<int> (pImpl->mXOffset/dx);
    pImpl->mHaveXLocation = true;
}

double Source3D::getLocationInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mX;
}

double Source3D::getOffsetInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mXOffset;
}

int Source3D::getCellInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mCellX;
}

bool Source3D::haveLocationInX() const noexcept
{
    return pImpl->mHaveXLocation;
}

/// y location
void Source3D::setLocationInY(const double y)
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    auto ny = pImpl->mGeometry.getNumberOfGridPointsInY();
    double dy = pImpl->mGeometry.getGridSpacingInY();
    double y0 = pImpl->mGeometry.getOriginInY();
    double y1 = y0 + static_cast<double>  (ny - 1)*dy;
    if (y < y0 || y > y1) 
    {
        throw std::invalid_argument("y = " + std::to_string(y)
                                  + " must be in range [" + std::to_string(y0)
                                  + "," + std::to_string(y1) + "]");
    }   
    pImpl->mY = y;
    pImpl->mYOffset = y - y0;
    pImpl->mCellY = static_cast<int> (pImpl->mYOffset/dy);
    pImpl->mHaveYLocation = true;
}

double Source3D::getLocationInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mY;
}

double Source3D::getOffsetInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mYOffset;
}

int Source3D::getCellInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mCellY;
}

bool Source3D::haveLocationInY() const noexcept
{
    return pImpl->mHaveYLocation;
}

/// z location
void Source3D::setLocationInZ(const double z)
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
    double dz = pImpl->mGeometry.getGridSpacingInZ();
    double z0 = pImpl->mGeometry.getOriginInZ();
    double z1 = z0 + static_cast<double>  (nz - 1)*dz;
    if (z < z0 || z > z1)
    {
        throw std::invalid_argument("z = " + std::to_string(z)
                                  + " must be in range [" + std::to_string(z0)
                                  + "," + std::to_string(z1) + "]");
    }
    pImpl->mZ = z;
    pImpl->mZOffset = z - z0;
    pImpl->mCellZ = static_cast<int> (pImpl->mZOffset/dz);
    pImpl->mHaveZLocation = true;
}

void Source3D::setZToFreeSurface()
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    double z0 = pImpl->mGeometry.getOriginInZ();
    setLocationInZ(z0);
}

double Source3D::getLocationInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZ;
}

double Source3D::getOffsetInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZOffset;
}

int Source3D::getCellInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCellZ;
}

bool Source3D::haveLocationInZ() const noexcept
{
    return pImpl->mHaveZLocation;
}

/// sdt::cout << source << std::endl;
std::ostream&
EikonalXX::operator<<(std::ostream &os, const Source3D &source)
{
    std::string result;
    if (source.haveLocationInX() &&
        source.haveLocationInY() &&
        source.haveLocationInZ())
    {
        result = "Source location: (x,y,z) = ("
               + std::to_string(source.getLocationInX()) + ","
               + std::to_string(source.getLocationInY()) + "," 
               + std::to_string(source.getLocationInZ()) + ")\n"
               + "Source offset: (x,y,z) = ("
               + std::to_string(source.getOffsetInX()) + ","
               + std::to_string(source.getOffsetInY()) + ","
               + std::to_string(source.getOffsetInZ()) + ")\n"
               + "Source cell: (iCellX,iCellY,iCellZ) = ("
               + std::to_string(source.getCellInX()) + ","
               + std::to_string(source.getCellInY()) + ","
               + std::to_string(source.getCellInZ()) + ")\n";
    }
    else
    {
        result = "Source location not completely specified";
    }
    return os << result;
}
