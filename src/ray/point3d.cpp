#include "eikonalxx/ray/point3d.hpp"

using namespace EikonalXX::Ray;

class Point3D::Point3DImpl
{
public:
    double mX{0};
    double mY{0};
    double mZ{0};
    bool mHaveX{false};
    bool mHaveY{false};
    bool mHaveZ{false};
};

/// Constructor
Point3D::Point3D() :
    pImpl(std::make_unique<Point3DImpl> ())
{
}

/// Constructor with (x,y,z)
Point3D::Point3D(const double x, const double y, const double z) :
    pImpl(std::make_unique<Point3DImpl> ()) 
{
    setPositionInX(x);
    setPositionInY(y);
    setPositionInZ(z);
}


/// Copy constructor
Point3D::Point3D(const Point3D &point)
{
    *this = point;
}

/// Move constructor
Point3D::Point3D(Point3D &&point) noexcept
{
    *this = std::move(point);
}

/// Copy assignment
Point3D& Point3D::operator=(const Point3D &point)
{
    if (&point == this){return *this;}
    pImpl = std::make_unique<Point3DImpl> (*point.pImpl);
    return *this;
}

/// Move assignment
Point3D& Point3D::operator=(Point3D &&point) noexcept
{
    if (&point == this){return *this;}
    pImpl = std::move(point.pImpl);
    return *this;
}

/// Position in x
void Point3D::setPositionInX(const double x) noexcept
{
    pImpl->mX = x;
    pImpl->mHaveX = true;
}

double Point3D::getPositionInX() const
{
    if (!havePositionInX()){throw std::runtime_error("x position not set");}
    return pImpl->mX;
}

bool Point3D::havePositionInX() const noexcept
{
    return pImpl->mHaveX;
}


/// Position in y
void Point3D::setPositionInY(const double y) noexcept
{
    pImpl->mY = y;
    pImpl->mHaveY = true;
}

double Point3D::getPositionInY() const
{
    if (!havePositionInY()){throw std::runtime_error("y position not set");}
    return pImpl->mY;
}

bool Point3D::havePositionInY() const noexcept
{
    return pImpl->mHaveY;
}

/// Position in z
void Point3D::setPositionInZ(const double z) noexcept
{
    pImpl->mZ = z;
    pImpl->mHaveZ = true;
}

double Point3D::getPositionInZ() const
{
    if (!havePositionInZ()){throw std::runtime_error("z position not set");}
    return pImpl->mZ;
}

bool Point3D::havePositionInZ() const noexcept
{
    return pImpl->mHaveZ;
}

/// Clear
void Point3D::clear() noexcept
{
    pImpl = std::make_unique<Point3DImpl> ();
}

/// Destructor
Point3D::~Point3D() = default;
