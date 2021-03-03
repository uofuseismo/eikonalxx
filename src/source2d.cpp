#include "eikonalxx/source2d.hpp"
#include "eikonalxx/geometry2d.hpp"

using namespace EikonalXX;

class Source2D::Source2DImpl
{
public:
    Geometry2D mGeometry;
    double mX = 0; 
    double mZ = 0;
    double mXOffset = 0;
    double mZOffset = 0;
    int mCellX = 0;
    int mCellZ = 0;
    bool mHaveXLocation = false;
    bool mHaveZLocation = false;
    bool mHaveGeometry = false;
};

/// Reset the class
void Source2D::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mX = 0;
    pImpl->mZ = 0;
    pImpl->mXOffset = 0;
    pImpl->mZOffset = 0;
    pImpl->mCellX = 0;
    pImpl->mCellZ = 0;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = false;
}

/// C'tor
Source2D::Source2D() :
    pImpl(std::make_unique<Source2DImpl> ())
{
}

/// Copy c'tor
Source2D::Source2D(const Source2D &source)
{
    *this = source;
}

/// Move c'tor
Source2D::Source2D(Source2D &&source) noexcept
{
    *this = std::move(source);
}

/// Copy assignment operator
Source2D& Source2D::operator=(const Source2D &source)
{
    if (&source == this){return *this;}
    pImpl = std::make_unique<Source2DImpl> (*source.pImpl);
    return *this;
}

/// Move assignment operator
Source2D& Source2D::operator=(Source2D &&source) noexcept
{
    if (&source == this){return *this;}
    pImpl = std::move(source.pImpl);
    return *this;
}

/// Destructor
Source2D::~Source2D() = default;

/// Geometry
void Source2D::setGeometry(const Geometry2D &geometry)
{
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x must be set");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("Grid spacing in z must be set");
    }
    if (!geometry.haveNumberOfGridPointsInX())
    {   
        throw std::invalid_argument("Number of grid points in x must be set");
    }
    if (!geometry.haveGridSpacingInZ())
    {   
        throw std::invalid_argument("Grid spacing in z must be set");
    }
    pImpl->mGeometry = geometry;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = true;
}

Geometry2D Source2D::getGeometry() const
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    return pImpl->mGeometry;
}

bool Source2D::haveGeometry() const noexcept
{
    return pImpl->mHaveGeometry;
}

/// x location
void Source2D::setLocationInX(const double x)
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

double Source2D::getLocationInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mX;
}

int Source2D::getSourceCellInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mCellX;
}

bool Source2D::haveLocationInX() const noexcept
{
    return pImpl->mHaveXLocation;
}

/// z location
void Source2D::setLocationInZ(const double z)
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

double Source2D::getLocationInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZ;
}

int Source2D::getSourceCellInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCellZ;
}

bool Source2D::haveLocationInZ() const noexcept
{
    return pImpl->mHaveZLocation;
}

