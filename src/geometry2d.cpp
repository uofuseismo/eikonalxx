#include <cmath>
#include <stdexcept>
#include <string>
#include "eikonalxx/geometry2d.hpp"

using namespace EikonalXX;

class Geometry2D::Geometry2DImpl
{
public: 
    double mDeltaX =-1;
    double mDeltaZ =-1;
    double mX0 = 0;
    double mZ0 = 0;
    int mGridX =-1;
    int mGridZ =-1;
};

/// C'tor
Geometry2D::Geometry2D() :
    pImpl(std::make_unique<Geometry2DImpl> ())
{
}

/// Copy c'tor
Geometry2D::Geometry2D(const Geometry2D &geometry)
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
    pImpl = std::make_unique<Geometry2DImpl> (*geometry.pImpl);
    return *this;
}

/// Move assingment
Geometry2D& Geometry2D::operator=(Geometry2D &&geometry) noexcept
{
    if (&geometry == this){return *this;}
    pImpl = std::move(geometry.pImpl);
    return *this;
}

/// Destructor
Geometry2D::~Geometry2D() = default;

/// Clears the class
void Geometry2D::clear() noexcept
{
    pImpl->mDeltaX =-1;
    pImpl->mDeltaZ =-1;
    pImpl->mX0 = 0;
    pImpl->mZ0 = 0;
    pImpl->mGridX =-1;
    pImpl->mGridZ =-1;
}

/// Set grid points in x
void Geometry2D::setNumberOfGridPointsInX(const int nx)
{
    if (nx < 3)
    {
        throw std::invalid_argument("Number of grid points in x = "
                                  + std::to_string(nx) + " must be at least 3");
    }
    pImpl->mGridX = nx;
}

/// Get grid points in x
int Geometry2D::getNumberOfGridPointsInX() const
{
    if (!haveNumberOfGridPointsInX())
    {
        throw std::runtime_error("Number of grid points in x not yet set");
    }
    return pImpl->mGridX;
}

/// Get cells in x
int Geometry2D::getNumberOfCellsInX() const
{
    return getNumberOfGridPointsInX() - 1;
}
 
/// Have grid points in x?
bool Geometry2D::haveNumberOfGridPointsInX() const noexcept
{
    return pImpl->mGridX > 0;
}

/// Set grid points in z
void Geometry2D::setNumberOfGridPointsInZ(const int nz)
{
    if (nz < 3)
    {
        throw std::invalid_argument("Number of grid points in z = "
                                  + std::to_string(nz) + " must be at least 3");
    }
    pImpl->mGridZ = nz;
}

/// Get grid points in z
int Geometry2D::getNumberOfGridPointsInZ() const
{
    if (!haveNumberOfGridPointsInZ())
    {
        throw std::runtime_error("Number of grid points in z not yet set");
    }
    return pImpl->mGridZ;
}

/// Get cells in z
int Geometry2D::getNumberOfCellsInZ() const
{
    return getNumberOfGridPointsInZ() - 1;
}

/// Have grid points in z?
bool Geometry2D::haveNumberOfGridPointsInZ() const noexcept
{
    return pImpl->mGridZ > 0;
}

/// Get number of grid points
int Geometry2D::getNumberOfGridPoints() const
{
    auto n = getNumberOfGridPointsInX()*getNumberOfGridPointsInZ();
    return n;
}

/// Get number of cells 
int Geometry2D::getNumberOfCells() const
{
    auto n = getNumberOfCellsInX()*getNumberOfCellsInZ();
    return n;
}

/// Set grid spacing in x
void Geometry2D::setGridSpacingInX(const double dx)
{
    if (dx <= 0)
    {
        throw std::invalid_argument("Grid spacing in x = " 
                                  + std::to_string(dx) + " must be positive");
    }
    pImpl->mDeltaX = dx;
}

/// Get grid spacing in x
double Geometry2D::getGridSpacingInX() const
{
    if (!haveGridSpacingInX())
    {
        throw std::runtime_error("Grid spacing in x not yet set");
    }
    return pImpl->mDeltaX;
}

/// Have grid spacing in x?
bool Geometry2D::haveGridSpacingInX() const noexcept
{
    return pImpl->mDeltaX > 0;
}

/// Grid spacing in z
void Geometry2D::setGridSpacingInZ(const double dz)
{
    if (dz <= 0)
    {
        throw std::invalid_argument("Grid spacing in z = "
                                  + std::to_string(dz) + " must be positive");
    }
    pImpl->mDeltaZ = dz;
}

/// Get grid spacing in z
double Geometry2D::getGridSpacingInZ() const
{
    if (!haveGridSpacingInZ())
    {
        throw std::runtime_error("Grid spacing in z not yet set");
    }
    return pImpl->mDeltaZ;
}

/// Have grid spacing in x?
bool Geometry2D::haveGridSpacingInZ() const noexcept
{
    return pImpl->mDeltaZ > 0;
}

/// X origin
void Geometry2D::setOriginInX(const double x0) noexcept
{
    pImpl->mX0 = x0;
}

double Geometry2D::getOriginInX() const noexcept
{
    return pImpl->mX0;
}

/// Z origin
void Geometry2D::setOriginInZ(const double z0) noexcept
{
    pImpl->mZ0 = z0;
}

double Geometry2D::getOriginInZ() const noexcept
{
    return pImpl->mZ0;
}

/// lhs == rhs
bool EikonalXX::operator==(const Geometry2D &lhs, const Geometry2D &rhs)
{
    // Preliminary check to avoid subsequent checks throwing an uninitialized 
    // error.
    if (lhs.haveNumberOfGridPointsInX() != rhs.haveNumberOfGridPointsInX())
    {
        return false;
    }
    if (lhs.haveNumberOfGridPointsInZ() != rhs.haveNumberOfGridPointsInZ())
    {
        return false;
    }
    if (lhs.haveGridSpacingInX() != rhs.haveGridSpacingInX())
    {
        return false;
    }
    if (lhs.haveGridSpacingInZ() != rhs.haveGridSpacingInZ())
    {
        return false;
    }
    // Check number of grid points
    if (lhs.haveNumberOfGridPointsInX())
    {
        if (lhs.getNumberOfGridPointsInX() != rhs.getNumberOfGridPointsInX())
        {
            return false;
        }
    }
    if (lhs.haveNumberOfGridPointsInZ())
    {
        if (lhs.getNumberOfGridPointsInZ() != rhs.getNumberOfGridPointsInZ())
        {
            return false;
        }
    }
    // Check the grid spacing
    if (lhs.haveGridSpacingInX())
    {
        if (std::abs(lhs.getGridSpacingInX() - rhs.getGridSpacingInX()) > 1.e-8)
        { 
            return false;
        }
    }
    if (lhs.haveGridSpacingInZ())
    {
        if (std::abs(lhs.getGridSpacingInZ() - rhs.getGridSpacingInZ()) > 1.e-8)
        {
            return false;
        }
    }
    // Check origin
    if (std::abs(lhs.getOriginInX() - rhs.getOriginInX()) > 1.e-7)
    {
        return false;
    }
    if (std::abs(lhs.getOriginInZ() - rhs.getOriginInZ()) > 1.e-7)
    {
        return false;
    }
    return true;
}

/// lhs != rhs
bool EikonalXX::operator!=(const Geometry2D &lhs, const Geometry2D &rhs)
{
    return !(lhs == rhs);
}
