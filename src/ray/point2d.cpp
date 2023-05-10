#include "eikonalxx/ray/point2d.hpp"

using namespace EikonalXX::Ray;

class Point2D::Point2DImpl
{
public:
    double mX{0};
    double mZ{0};
    bool mHaveX{false};
    bool mHaveZ{false};
};

/// Constructor
Point2D::Point2D() :
    pImpl(std::make_unique<Point2DImpl> ())
{
}

/// Copy constructor
Point2D::Point2D(const Point2D &point)
{
    *this = point;
}

/// Move constructor
Point2D::Point2D(Point2D &&point) noexcept
{
    *this = std::move(point);
}

/// Copy assignment
Point2D& Point2D::operator=(const Point2D &point)
{
    if (&point == this){return *this;}
    pImpl = std::make_unique<Point2DImpl> (*point.pImpl);
    return *this;
}

/// Move assignment
Point2D& Point2D::operator=(Point2D &&point) noexcept
{
    if (&point == this){return *this;}
    pImpl = std::move(point.pImpl);
    return *this;
}

/// Position in x
void Point2D::setPositionInX(const double x) noexcept
{
    pImpl->mX = x;
    pImpl->mHaveX = true;
}

double Point2D::getPositionInX() const
{
    if (!havePositionInX()){throw std::runtime_error("x position not set");}
    return pImpl->mX;
}

bool Point2D::havePositionInX() const noexcept
{
    return pImpl->mHaveX;
}

/// Position in z
void Point2D::setPositionInZ(const double z) noexcept
{
    pImpl->mZ = z;
    pImpl->mHaveZ = true;
}

double Point2D::getPositionInZ() const
{
    if (!havePositionInZ()){throw std::runtime_error("z position not set");}
    return pImpl->mZ;
}

bool Point2D::havePositionInZ() const noexcept
{
    return pImpl->mHaveZ;
}

/// Clear
void Point2D::clear() noexcept
{
    pImpl = std::make_unique<Point2DImpl> ();
}

/// Destructor
Point2D::~Point2D() = default;
