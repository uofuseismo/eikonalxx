#include "eikonalxx/station3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX;

class Station3D::Station3DImpl
{
public:
    void updateCell()
    {
        if (mCell < 0)
        {
            if (mHaveXLocation && mHaveYLocation && mHaveZLocation)
            {
                auto nCellX = mGeometry.getNumberOfCellsInX();
                auto nCellY = mGeometry.getNumberOfCellsInY();
                mCell = gridToIndex(nCellX, nCellY, mCellX, mCellY, mCellZ);
            }
        }
    }
    Geometry3D mGeometry;
    std::string mName;
    double mX = 0; 
    double mY = 0;
    double mZ = 0;
    double mXOffset = 0;
    double mYOffset = 0;
    double mZOffset = 0;
    int mCellX = 0;
    int mCellY = 0;
    int mCellZ = 0;
    int mCell =-1;
    bool mHaveXLocation = false;
    bool mHaveYLocation = false;
    bool mHaveZLocation = false;
    bool mHaveGeometry = false;
};

/// Reset the class
void Station3D::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mName.clear();
    pImpl->mX = 0;
    pImpl->mY = 0;
    pImpl->mZ = 0;
    pImpl->mXOffset = 0;
    pImpl->mYOffset = 0;
    pImpl->mZOffset = 0;
    pImpl->mCellX = 0;
    pImpl->mCellY = 0;
    pImpl->mCellZ = 0;
    pImpl->mCell =-1;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveYLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = false;
}

/// C'tor
Station3D::Station3D() :
    pImpl(std::make_unique<Station3DImpl> ())
{
}

/// Copy c'tor
Station3D::Station3D(const Station3D &station)
{
    *this = station;
}

/// Move c'tor
Station3D::Station3D(Station3D &&station) noexcept
{
    *this = std::move(station);
}

/// Copy assignment operator
Station3D& Station3D::operator=(const Station3D &station)
{
    if (&station == this){return *this;}
    pImpl = std::make_unique<Station3DImpl> (*station.pImpl);
    return *this;
}

/// Move assignment operator
Station3D& Station3D::operator=(Station3D &&station) noexcept
{
    if (&station == this){return *this;}
    pImpl = std::move(station.pImpl);
    return *this;
}

/// Destructor
Station3D::~Station3D() = default;

/// Geometry
void Station3D::setGeometry(const Geometry3D &geometry)
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

Geometry3D Station3D::getGeometry() const
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    return pImpl->mGeometry;
}

bool Station3D::haveGeometry() const noexcept
{
    return pImpl->mHaveGeometry;
}

/// Set the station name
void Station3D::setName(const std::string &name) noexcept
{
    pImpl->mName = name;
}

/// Get the station name
std::string Station3D::getName() const noexcept
{
    return pImpl->mName;
}

/// Cell
int Station3D::getCell() const
{
    if (pImpl->mCell < 0)
    {
        if (!haveLocationInX())
        {
            throw std::runtime_error("x location not yet set");
        }
        if (!haveLocationInY())
        {
            throw std::runtime_error("y location not yet set");
        }
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCell;
}

/// x location
void Station3D::setLocationInX(const double x)
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
    pImpl->mCell =-1;
    pImpl->mX = x;
    pImpl->mXOffset = x - x0;
    pImpl->mCellX = static_cast<int> (pImpl->mXOffset/dx);
    pImpl->mHaveXLocation = true;
    pImpl->updateCell();
}

double Station3D::getLocationInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mX;
}

double Station3D::getOffsetInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mXOffset;
}

int Station3D::getCellInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mCellX;
}

bool Station3D::haveLocationInX() const noexcept
{
    return pImpl->mHaveXLocation;
}

/// y location
void Station3D::setLocationInY(const double y)
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
    pImpl->mCell =-1;
    pImpl->mY = y;
    pImpl->mYOffset = y - y0;
    pImpl->mCellY = static_cast<int> (pImpl->mYOffset/dy);
    pImpl->mHaveYLocation = true;
    pImpl->updateCell();
}

double Station3D::getLocationInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mY;
}

double Station3D::getOffsetInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mYOffset;
}

int Station3D::getCellInY() const
{
    if (!haveLocationInY())
    {
        throw std::runtime_error("y location not yet set");
    }
    return pImpl->mCellY;
}

bool Station3D::haveLocationInY() const noexcept
{
    return pImpl->mHaveYLocation;
}

/// z location
void Station3D::setLocationInZ(const double z)
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
    pImpl->mCell =-1;
    pImpl->mZ = z;
    pImpl->mZOffset = z - z0;
    pImpl->mCellZ = static_cast<int> (pImpl->mZOffset/dz);
    pImpl->mHaveZLocation = true;
    pImpl->updateCell();
}

void Station3D::setZToFreeSurface()
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    double z0 = pImpl->mGeometry.getOriginInZ();
    setLocationInZ(z0);
}

double Station3D::getLocationInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZ;
}

double Station3D::getOffsetInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZOffset;
}

int Station3D::getCellInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCellZ;
}

bool Station3D::haveLocationInZ() const noexcept
{
    return pImpl->mHaveZLocation;
}

/// sdt::cout << station << std::endl;
std::ostream&
EikonalXX::operator<<(std::ostream &os, const Station3D &station)
{
    std::string result;
    if (station.haveLocationInX() &&
        station.haveLocationInY() &&
        station.haveLocationInZ())
    {
        result = "Station " + station.getName() + " location: (x,y,z) = ("
               + std::to_string(station.getLocationInX()) + ","
               + std::to_string(station.getLocationInY()) + "," 
               + std::to_string(station.getLocationInZ()) + ")\n"
               + "Station offset: (x,y,z) = ("
               + std::to_string(station.getOffsetInX()) + ","
               + std::to_string(station.getOffsetInY()) + ","
               + std::to_string(station.getOffsetInZ()) + ")\n"
               + "Station cell: (iCellX,iCellY,iCellZ) = ("
               + std::to_string(station.getCellInX()) + ","
               + std::to_string(station.getCellInY()) + ","
               + std::to_string(station.getCellInZ()) + ")\n";
    }
    else
    {
        result = "Station location not completely specified";
    }
    return os << result;
}
