#include "eikonalxx/station2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX;

class Station2D::Station2DImpl
{
public:
    void updateCell()
    {
        if (mCell < 0)
        {
            if (mHaveXLocation && mHaveZLocation)
            {
                auto nCellX = mGeometry.getNumberOfCellsInX();
                mCell = ::gridToIndex(nCellX, mCellX, mCellZ);
            }
        }
    }
    Geometry2D mGeometry;
    std::string mName;
    double mX{0};
    double mZ{0};
    double mXOffset{0};
    double mZOffset{0};
    int mCellX{0};
    int mCellZ{0};
    int mCell{-1};
    bool mHaveXLocation{false};
    bool mHaveZLocation{false};
    bool mHaveGeometry{false};
};

/// Reset the class
void Station2D::clear() noexcept
{
    pImpl = std::make_unique<Station2DImpl> ();
}

/// C'tor
Station2D::Station2D() :
    pImpl(std::make_unique<Station2DImpl> ())
{
}

/// Copy c'tor
Station2D::Station2D(const Station2D &station)
{
    *this = station;
}

/// Move c'tor
Station2D::Station2D(Station2D &&station) noexcept
{
    *this = std::move(station);
}

/// Copy assignment operator
Station2D& Station2D::operator=(const Station2D &station)
{
    if (&station == this){return *this;}
    pImpl = std::make_unique<Station2DImpl> (*station.pImpl);
    return *this;
}

/// Move assignment operator
Station2D& Station2D::operator=(Station2D &&station) noexcept
{
    if (&station == this){return *this;}
    pImpl = std::move(station.pImpl);
    return *this;
}

/// Destructor
Station2D::~Station2D() = default;

/// Geometry
void Station2D::setGeometry(const Geometry2D &geometry)
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
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Number of grid points in z must be set");
    }
    pImpl->mGeometry = geometry;
    pImpl->mHaveXLocation = false;
    pImpl->mHaveZLocation = false;
    pImpl->mHaveGeometry = true;
}

Geometry2D Station2D::getGeometry() const
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    return pImpl->mGeometry;
}

bool Station2D::haveGeometry() const noexcept
{
    return pImpl->mHaveGeometry;
}

/// Set the station name
void Station2D::setName(const std::string &name) noexcept
{
    pImpl->mName = name;
}

/// Get the station name
std::string Station2D::getName() const noexcept
{
    return pImpl->mName;
}

/// Cell
int Station2D::getCell() const
{
    if (pImpl->mCell < 0)
    {
        if (!haveLocationInX())
        {
            throw std::runtime_error("x location not yet set");
        }
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCell;
}

/// x location
void Station2D::setLocationInX(const double x)
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    auto nCellX = pImpl->mGeometry.getNumberOfCellsInX();
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
    pImpl->mCellX = std::min(nCellX - 1, static_cast<int> (pImpl->mXOffset/dx));
    pImpl->mHaveXLocation = true;
    pImpl->updateCell();
}

double Station2D::getLocationInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mX;
}

double Station2D::getOffsetInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mXOffset;
}

int Station2D::getCellInX() const
{
    if (!haveLocationInX())
    {
        throw std::runtime_error("x location not yet set");
    }
    return pImpl->mCellX;
}

bool Station2D::haveLocationInX() const noexcept
{
    return pImpl->mHaveXLocation;
}

/// z location
void Station2D::setLocationInZ(const double z)
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
    auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
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
    pImpl->mCellZ = std::min(nCellZ - 1, static_cast<int> (pImpl->mZOffset/dz));
    pImpl->mHaveZLocation = true;
    pImpl->updateCell();
}

void Station2D::setZToFreeSurface()
{
    if (!haveGeometry()){throw std::runtime_error("Geometry not yet set");}
    double z0 = pImpl->mGeometry.getOriginInZ();
    setLocationInZ(z0);
}

double Station2D::getLocationInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZ;
}

double Station2D::getOffsetInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mZOffset;
}

int Station2D::getCellInZ() const
{
    if (!haveLocationInZ())
    {
        throw std::runtime_error("z location not yet set");
    }
    return pImpl->mCellZ;
}

bool Station2D::haveLocationInZ() const noexcept
{
    return pImpl->mHaveZLocation;
}

/// std::cout << station << std::endl;
std::ostream&
EikonalXX::operator<<(std::ostream &os, const Station2D &station)
{
    std::string result;
    if (station.haveLocationInX() && station.haveLocationInZ())
    {
        result = "Station " + station.getName() + " location: (x,z) = ("
               + std::to_string(station.getLocationInX()) + ","
               + std::to_string(station.getLocationInZ()) + ")\n"
               + "Station offset: (x,z) = ("
               + std::to_string(station.getOffsetInX()) + ","
               + std::to_string(station.getOffsetInZ()) + ")\n"
               + "Station cell: (iCellX,iCellY,iCellZ) = ("
               + std::to_string(station.getCellInX()) + ","
               + std::to_string(station.getCellInZ()) + ")\n";
    }
    else
    {
        result = "Station location not completely specified";
    }
    return os << result;
}
