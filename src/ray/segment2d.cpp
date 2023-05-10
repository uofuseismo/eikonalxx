#include <cmath>
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"

using namespace EikonalXX::Ray;

namespace
{

double computeLength(const Point2D &a, const Point2D &b)
{
    auto dx = b.getPositionInX() - a.getPositionInX();
    auto dz = b.getPositionInZ() - a.getPositionInZ();
    return std::hypot(dx, dz);
}

}

class Segment2D::Segment2DImpl
{
public:
    Point2D mStartPoint;
    Point2D mEndPoint;
    double mLength{0};
    double mVelocity{0};
    int mVelocityModelCellIndex{-1};
    bool mHaveStartAndEndPoint{false};
};

/// Constructor
Segment2D::Segment2D() :
    pImpl(std::make_unique<Segment2DImpl> ())
{
}

/// Copy constructor
Segment2D::Segment2D(const Segment2D &segment)
{
    *this = segment;
}

/// Move constructor
Segment2D::Segment2D(Segment2D &&segment) noexcept
{
    *this = std::move(segment);
}

/// Copy assignment
Segment2D& Segment2D::operator=(const Segment2D &segment)
{
    if (&segment == this){return *this;}
    pImpl = std::make_unique<Segment2DImpl> (*segment.pImpl);
    return *this;
}

/// Move assignment
Segment2D& Segment2D::operator=(Segment2D &&segment) noexcept
{
    if (&segment == this){return *this;}
    pImpl = std::move(segment.pImpl);
    return *this;
}

/// Reset class
void Segment2D::clear() noexcept
{
    pImpl = std::make_unique<Segment2DImpl> ();
}

/// Destructor
Segment2D::~Segment2D() = default;

/// Start and end point
void Segment2D::setStartAndEndPoint(
    const std::pair<Point2D, Point2D> &startAndEndPoint)
{
    if (!startAndEndPoint.first.havePositionInX())
    {
        throw std::invalid_argument("Start x position not set");
    }
    if (!startAndEndPoint.first.havePositionInZ())
    {
        throw std::invalid_argument("Start z position not set");
    }
    if (!startAndEndPoint.second.havePositionInX())
    {
        throw std::invalid_argument("End x position not set");
    }
    if (!startAndEndPoint.second.havePositionInZ())
    {
        throw std::invalid_argument("End z position not set");
    }
    pImpl->mStartPoint = startAndEndPoint.first;
    pImpl->mEndPoint = startAndEndPoint.second;
    pImpl->mLength = ::computeLength(pImpl->mStartPoint, pImpl->mEndPoint);
    pImpl->mHaveStartAndEndPoint = true;
}

Point2D Segment2D::getStartPoint() const
{
    if (!haveStartAndEndPoint())
    {   
        throw std::runtime_error("Start point not set");
    }   
    return pImpl->mStartPoint;
}

Point2D Segment2D::getEndPoint() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("End point not set");
    }
    return pImpl->mEndPoint;
}

double Segment2D::getLength() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("Start and end point not set");
    }
    return pImpl->mLength;
}

bool Segment2D::haveStartAndEndPoint() const noexcept
{
    return pImpl->mHaveStartAndEndPoint;
}

/// Velocity
void Segment2D::setVelocity(const double velocity)
{
    if (velocity <= 0){throw std::runtime_error("Velocity must be positive");}
    pImpl->mVelocity = velocity;
}

double Segment2D::getVelocity() const
{
    if (!haveVelocity()){throw std::runtime_error("Velocity not set");}
    return pImpl->mVelocity;
}

double Segment2D::getSlowness() const
{
    return 1./getVelocity();
}

bool Segment2D::haveVelocity() const noexcept
{
    return (pImpl->mVelocity > 0);
}

/// Travel time
double Segment2D::getTravelTime() const
{
    return getSlowness()*getLength();
}

/// Cell index
void Segment2D::setVelocityModelCellIndex(const int cellIndex)
{
    if (cellIndex < 0)
    {
        throw std::invalid_argument("Cell index must be positive");
    }
    pImpl->mVelocityModelCellIndex = cellIndex;
}

int Segment2D::getVelocityModelCellIndex() const
{
    if (!haveVelocityModelCellIndex())
    {
        throw std::runtime_error("Velocity model cell index not set");
    }
    return pImpl->mVelocityModelCellIndex;
}

bool Segment2D::haveVelocityModelCellIndex() const noexcept
{
    return (pImpl->mVelocityModelCellIndex >= 0);
}
