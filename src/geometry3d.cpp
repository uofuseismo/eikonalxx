#include <cmath>
#include <stdexcept>
#include <string>
#include "eikonalxx/geometry3d.hpp"

using namespace EikonalXX;

class Geometry3D::Geometry3DImpl
{
public: 
    double mDeltaX =-1;
    double mDeltaY =-1;
    double mDeltaZ =-1;
    double mX0 = 0;
    double mY0 = 0;
    double mZ0 = 0;
    int mGridX =-1;
    int mGridY =-1;
    int mGridZ =-1;
};

/// C'tor
Geometry3D::Geometry3D() :
    pImpl(std::make_unique<Geometry3DImpl> ())
{
}

/// Copy c'tor
Geometry3D::Geometry3D(const Geometry3D &geometry)
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
    pImpl = std::make_unique<Geometry3DImpl> (*geometry.pImpl);
    return *this;
}

/// Move assingment
Geometry3D& Geometry3D::operator=(Geometry3D &&geometry) noexcept
{
    if (&geometry == this){return *this;}
    pImpl = std::move(geometry.pImpl);
    return *this;
}

/// Destructor
Geometry3D::~Geometry3D() = default;

/// Clears the class
void Geometry3D::clear() noexcept
{
    pImpl->mDeltaX =-1;
    pImpl->mDeltaY =-1;
    pImpl->mDeltaZ =-1;
    pImpl->mX0 = 0;
    pImpl->mY0 = 0;
    pImpl->mZ0 = 0;
    pImpl->mGridX =-1;
    pImpl->mGridY =-1;
    pImpl->mGridZ =-1;
}

/// Set grid points in x
void Geometry3D::setNumberOfGridPointsInX(const int nx)
{
    if (nx < 3)
    {
        throw std::invalid_argument("Number of grid points in x = "
                                  + std::to_string(nx) + " must be at least 3");
    }
    pImpl->mGridX = nx;
}

/// Get grid points in x
int Geometry3D::getNumberOfGridPointsInX() const
{
    if (!haveNumberOfGridPointsInX())
    {
        throw std::runtime_error("Number of grid points in x not yet set");
    }
    return pImpl->mGridX;
}

/// Get cells in x
int Geometry3D::getNumberOfCellsInX() const
{
    return getNumberOfGridPointsInX() - 1;
}
 
/// Have grid points in x?
bool Geometry3D::haveNumberOfGridPointsInX() const noexcept
{
    return pImpl->mGridX > 0;
}

/// Set grid points in y
void Geometry3D::setNumberOfGridPointsInY(const int ny)
{
    if (ny < 3)
    {
        throw std::invalid_argument("Number of grid points in y = "
                                  + std::to_string(ny) + " must be at least 3");
    }
    pImpl->mGridY = ny;
}

/// Get grid points in y
int Geometry3D::getNumberOfGridPointsInY() const
{
    if (!haveNumberOfGridPointsInY())
    {
        throw std::runtime_error("Number of grid points in y not yet set");
    }
    return pImpl->mGridY;
}

/// Get cells in y
int Geometry3D::getNumberOfCellsInY() const
{
    return getNumberOfGridPointsInY() - 1;
}

/// Have grid points in y?
bool Geometry3D::haveNumberOfGridPointsInY() const noexcept
{
    return pImpl->mGridY > 0;
}


/// Set grid points in z
void Geometry3D::setNumberOfGridPointsInZ(const int nz)
{
    if (nz < 3)
    {
        throw std::invalid_argument("Number of grid points in z = "
                                  + std::to_string(nz) + " must be at least 3");
    }
    pImpl->mGridZ = nz;
}

/// Get grid points in z
int Geometry3D::getNumberOfGridPointsInZ() const
{
    if (!haveNumberOfGridPointsInZ())
    {
        throw std::runtime_error("Number of grid points in z not yet set");
    }
    return pImpl->mGridZ;
}

/// Get cells in z
int Geometry3D::getNumberOfCellsInZ() const
{
    return getNumberOfGridPointsInZ() - 1;
}

/// Have grid points in z?
bool Geometry3D::haveNumberOfGridPointsInZ() const noexcept
{
    return pImpl->mGridZ > 0;
}

/// Get number of grid points
int Geometry3D::getNumberOfGridPoints() const
{
    auto n = getNumberOfGridPointsInX()
            *getNumberOfGridPointsInY()
            *getNumberOfGridPointsInZ();
    return n;
}

/// Get number of cells 
int Geometry3D::getNumberOfCells() const
{
    auto n = getNumberOfCellsInX()
            *getNumberOfCellsInY()
            *getNumberOfCellsInZ();
    return n;
}

/// Set grid spacing in x
void Geometry3D::setGridSpacingInX(const double dx)
{
    if (dx <= 0)
    {
        throw std::invalid_argument("Grid spacing in x = " 
                                  + std::to_string(dx) + " must be positive");
    }
    pImpl->mDeltaX = dx;
}

/// Get grid spacing in x
double Geometry3D::getGridSpacingInX() const
{
    if (!haveGridSpacingInX())
    {
        throw std::runtime_error("Grid spacing in x not yet set");
    }
    return pImpl->mDeltaX;
}

/// Have grid spacing in x?
bool Geometry3D::haveGridSpacingInX() const noexcept
{
    return pImpl->mDeltaX > 0;
}

/// Grid spacing in y
void Geometry3D::setGridSpacingInY(const double dy)
{
    if (dy <= 0)
    {
        throw std::invalid_argument("Grid spacing in y = "
                                  + std::to_string(dy) + " must be positive");
    }
    pImpl->mDeltaY = dy;
}

/// Get grid spacing in y
double Geometry3D::getGridSpacingInY() const
{
    if (!haveGridSpacingInY())
    {   
        throw std::runtime_error("Grid spacing in y not yet set");
    }   
    return pImpl->mDeltaY;
}

/// Have grid spacing in y?
bool Geometry3D::haveGridSpacingInY() const noexcept
{
    return pImpl->mDeltaY > 0;
}

/// Grid spacing in z
void Geometry3D::setGridSpacingInZ(const double dz)
{
    if (dz <= 0)
    {
        throw std::invalid_argument("Grid spacing in z = "
                                  + std::to_string(dz) + " must be positive");
    }
    pImpl->mDeltaZ = dz;
}

/// Get grid spacing in z
double Geometry3D::getGridSpacingInZ() const
{
    if (!haveGridSpacingInZ())
    {
        throw std::runtime_error("Grid spacing in z not yet set");
    }
    return pImpl->mDeltaZ;
}

/// Have grid spacing in x?
bool Geometry3D::haveGridSpacingInZ() const noexcept
{
    return pImpl->mDeltaZ > 0;
}

/// X origin
void Geometry3D::setOriginInX(const double x0) noexcept
{
    pImpl->mX0 = x0;
}

double Geometry3D::getOriginInX() const noexcept
{
    return pImpl->mX0;
}

/// Y origin
void Geometry3D::setOriginInY(const double y0) noexcept
{
    pImpl->mY0 = y0;
}

double Geometry3D::getOriginInY() const noexcept
{
    return pImpl->mY0;
}

/// Z origin
void Geometry3D::setOriginInZ(const double z0) noexcept
{
    pImpl->mZ0 = z0;
}

double Geometry3D::getOriginInZ() const noexcept
{
    return pImpl->mZ0;
}

/// lhs == rhs
bool EikonalXX::operator==(const Geometry3D &lhs, const Geometry3D &rhs)
{
    // Preliminary check to avoid subsequent checks throwing an uninitialized 
    // error.
    if (lhs.haveNumberOfGridPointsInX() != rhs.haveNumberOfGridPointsInX())
    {
        return false;
    }
    if (lhs.haveNumberOfGridPointsInY() != rhs.haveNumberOfGridPointsInY())
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
    if (lhs.haveGridSpacingInY() != rhs.haveGridSpacingInY())
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
    if (lhs.haveNumberOfGridPointsInY())
    {   
        if (lhs.getNumberOfGridPointsInY() != rhs.getNumberOfGridPointsInY())
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
    if (lhs.haveGridSpacingInY())
    {   
        if (std::abs(lhs.getGridSpacingInY() - rhs.getGridSpacingInY()) > 1.e-8)
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
    if (std::abs(lhs.getOriginInY() - rhs.getOriginInY()) > 1.e-7)
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
bool EikonalXX::operator!=(const Geometry3D &lhs, const Geometry3D &rhs)
{
    return !(lhs == rhs);
}
