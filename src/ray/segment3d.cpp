#include <cmath>
#include "eikonalxx/ray/segment3d.hpp"
#include "eikonalxx/ray/point3d.hpp"

using namespace EikonalXX::Ray;

namespace
{

double computeLength(const Point3D &a, const Point3D &b)
{
    auto dx = b.getPositionInX() - a.getPositionInX();
    auto dy = b.getPositionInY() - a.getPositionInY();
    auto dz = b.getPositionInZ() - a.getPositionInZ();
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

}

class Segment3D::Segment3DImpl
{
public:
    Point3D mStartPoint;
    Point3D mEndPoint;
    double mLength{0};
    double mVelocity{0};
    int mVelocityModelCellIndex{-1};
    bool mHaveStartAndEndPoint{false};
};

/// Constructor
Segment3D::Segment3D() :
    pImpl(std::make_unique<Segment3DImpl> ())
{
}

/// Copy constructor
Segment3D::Segment3D(const Segment3D &segment)
{
    *this = segment;
}

/// Move constructor
Segment3D::Segment3D(Segment3D &&segment) noexcept
{
    *this = std::move(segment);
}

/// Copy assignment
Segment3D& Segment3D::operator=(const Segment3D &segment)
{
    if (&segment == this){return *this;}
    pImpl = std::make_unique<Segment3DImpl> (*segment.pImpl);
    return *this;
}

/// Move assignment
Segment3D& Segment3D::operator=(Segment3D &&segment) noexcept
{
    if (&segment == this){return *this;}
    pImpl = std::move(segment.pImpl);
    return *this;
}

/// Reset class
void Segment3D::clear() noexcept
{
    pImpl = std::make_unique<Segment3DImpl> ();
}

/// Destructor
Segment3D::~Segment3D() = default;

/// Start and end point
void Segment3D::setStartAndEndPoint(
    const std::pair<Point3D, Point3D> &startAndEndPoint)
{
    if (!startAndEndPoint.first.havePositionInX())
    {
        throw std::invalid_argument("Start x position not set");
    }
    if (!startAndEndPoint.first.havePositionInY())
    {
        throw std::invalid_argument("Start y position not set");
    }
    if (!startAndEndPoint.first.havePositionInZ())
    {
        throw std::invalid_argument("Start z position not set");
    }
    if (!startAndEndPoint.second.havePositionInX())
    {
        throw std::invalid_argument("End x position not set");
    }
    if (!startAndEndPoint.second.havePositionInY())
    {
        throw std::invalid_argument("End y position not set");
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

Point3D Segment3D::getStartPoint() const
{
    if (!haveStartAndEndPoint())
    {   
        throw std::runtime_error("Start point not set");
    }   
    return pImpl->mStartPoint;
}

Point3D Segment3D::getEndPoint() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("End point not set");
    }
    return pImpl->mEndPoint;
}

double Segment3D::getLength() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("Start and end point not set");
    }
    return pImpl->mLength;
}

bool Segment3D::haveStartAndEndPoint() const noexcept
{
    return pImpl->mHaveStartAndEndPoint;
}

/// Velocity
void Segment3D::setVelocity(const double velocity)
{
    if (velocity <= 0){throw std::runtime_error("Velocity must be positive");}
    pImpl->mVelocity = velocity;
}

double Segment3D::getVelocity() const
{
    if (!haveVelocity()){throw std::runtime_error("Velocity not set");}
    return pImpl->mVelocity;
}

double Segment3D::getSlowness() const
{
    return 1./getVelocity();
}

bool Segment3D::haveVelocity() const noexcept
{
    return (pImpl->mVelocity > 0);
}

/// Travel time
double Segment3D::getTravelTime() const
{
    return getSlowness()*getLength();
}

/// Cell index
void Segment3D::setVelocityModelCellIndex(const int cellIndex)
{
    if (cellIndex < 0)
    {
        throw std::invalid_argument("Cell index must be positive");
    }
    pImpl->mVelocityModelCellIndex = cellIndex;
}

int Segment3D::getVelocityModelCellIndex() const
{
    if (!haveVelocityModelCellIndex())
    {
        throw std::runtime_error("Velocity model cell index not set");
    }
    return pImpl->mVelocityModelCellIndex;
}

bool Segment3D::haveVelocityModelCellIndex() const noexcept
{
    return (pImpl->mVelocityModelCellIndex >= 0);
}
