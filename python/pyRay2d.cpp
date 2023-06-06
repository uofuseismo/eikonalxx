#include "include/pyRay2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"

using namespace PEikonalXX::Ray;
    
Point2D::Point2D() :
    pImpl(std::make_unique<EikonalXX::Ray::Point2D> ())
{
}

Point2D::Point2D(const Point2D &point)
{
    *this = point;
}

Point2D::Point2D(const EikonalXX::Ray::Point2D &point)
{
    *this = point;
}

Point2D::Point2D(Point2D &&point) noexcept
{
    *this = std::move(point);
}

Point2D& Point2D::operator=(const Point2D &point)
{
    if (&point == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Ray::Point2D> (*point.pImpl);
    return *this;
}

Point2D& Point2D::operator=(const EikonalXX::Ray::Point2D &point)
{
    pImpl = std::make_unique<EikonalXX::Ray::Point2D> (point);
    return *this;
}

Point2D& Point2D::operator=(Point2D &&point) noexcept
{
    if (&point == this){return *this;}
    pImpl = std::move(point.pImpl);
    return *this;
}
 
void Point2D::setPositionInX(double x) noexcept
{
    pImpl->setPositionInX(x);
}

double Point2D::getPositionInX() const
{
    return pImpl->getPositionInX();
}

bool Point2D::havePositionInX() const noexcept
{
    return pImpl->havePositionInX();
}

void Point2D::setPositionInZ(double z) noexcept
{
    pImpl->setPositionInZ(z);
}

double Point2D::getPositionInZ() const
{
    return pImpl->getPositionInZ();
}

bool Point2D::havePositionInZ() const noexcept
{
    return pImpl->havePositionInZ();
}

EikonalXX::Ray::Point2D& Point2D::getNativeClassReference() const
{
    return *pImpl;
}

Point2D::~Point2D() = default;

void Point2D::clear() noexcept
{
    pImpl->clear();
}

///--------------------------------------------------------------------------///
Segment2D::Segment2D() :
    pImpl(std::make_unique<EikonalXX::Ray::Segment2D> ())
{
}

Segment2D::Segment2D(const Segment2D &segment)
{
    *this = segment;
}

Segment2D::Segment2D(const EikonalXX::Ray::Segment2D &segment)
{
    *this = segment;
}

Segment2D::Segment2D(Segment2D &&segment) noexcept
{
    *this = std::move(segment);
}

Segment2D& Segment2D::operator=(const Segment2D &segment)
{
    if (&segment == this){return *this;}
    pImpl = std::make_unique<EikonalXX::Ray::Segment2D> (*segment.pImpl);
    return *this;
}

Segment2D& Segment2D::operator=(const EikonalXX::Ray::Segment2D &segment)
{
    pImpl = std::make_unique<EikonalXX::Ray::Segment2D> (segment);
    return *this;
}

Segment2D& Segment2D::operator=(Segment2D &&segment) noexcept
{
    if (&segment == this){return *this;}
    pImpl = std::move(segment.pImpl);
    return *this;
}

EikonalXX::Ray::Segment2D& Segment2D::getNativeClassReference() const
{
    return *pImpl;
}


void Segment2D::setStartAndEndPoint(const std::pair<Point2D, Point2D> &startAndEndPoint)
{
    auto startPoint = startAndEndPoint.first.getNativeClassReference();
    auto endPoint   = startAndEndPoint.second.getNativeClassReference();
    pImpl->setStartAndEndPoint(std::pair {startPoint, endPoint});
}

Point2D Segment2D::getStartPoint() const
{
    Point2D result{pImpl->getStartPoint()};
    return result;
}

Point2D Segment2D::getEndPoint() const
{
    Point2D result{pImpl->getEndPoint()};
    return result;
}

bool Segment2D::haveStartAndEndPoint() const noexcept
{
    return pImpl->haveStartAndEndPoint();
}

double Segment2D::getLength() const
{
    return pImpl->getLength();
}

void Segment2D::setVelocity(double velocity)
{
    pImpl->setVelocity(velocity);
}

void Segment2D::setSlowness(double slowness)
{
    pImpl->setSlowness(slowness);
}

double Segment2D::getVelocity() const
{
    return pImpl->getVelocity();
}

double Segment2D::getSlowness() const
{
    return pImpl->getSlowness();
}

bool Segment2D::haveVelocity() const noexcept
{
    return pImpl->haveVelocity();
}

double Segment2D::getTravelTime() const
{
    return pImpl->getTravelTime();
}

void Segment2D::setVelocityModelCellIndex(int index)
{
    pImpl->setVelocityModelCellIndex(index);
}

int Segment2D::getVelocityModelCellIndex() const
{
    return pImpl->getVelocityModelCellIndex();
}

bool Segment2D::haveVelocityModelCellIndex() const noexcept
{
    return pImpl->haveVelocityModelCellIndex();
}

Segment2D::~Segment2D() = default;

void Segment2D::clear() noexcept
{
    pImpl->clear();
}

///--------------------------------------------------------------------------///

/// Initialize module 
void PEikonalXX::Ray::initializeRay2D(pybind11::module &m)
{
    pybind11::class_<PEikonalXX::Ray::Point2D> point(m, "Point2D");
    point.def(pybind11::init<> ());
    point.doc() = R"""(
This defines a 2D Cartesian point along a ray.

Properties
----------
x : float
    The point's position (meters) in x 
z : float
    The point's position (meters) in z.
)""";
    point.def("__copy__", [](const Point2D &self)
    {
        return Point2D(self);
    });
    point.def_property("x",
                       &Point2D::getPositionInX,
                       &Point2D::setPositionInX);
    point.def_property("z",
                       &Point2D::getPositionInZ,
                       &Point2D::setPositionInZ);
    point.def("clear",
              &Point2D::clear,
              "Resets the class.");
    pybind11::class_<PEikonalXX::Ray::Segment2D> segment(m, "Segment2D");
    segment.def(pybind11::init<> ());
    segment.doc() = R"""(
This defines a 2D ray segment which is defined by a start and end point
as well as a velocity.

Properties
----------
velocity : float
    Sets the seismic velocity in m/s through which the ray segment passes.
velocity_model_cell_index : int
    The velocity model's cell index through which this ray passes.

Read-Only Properties
--------------------
start_point : Point2D
    The start point of the ray segment.
end_point : Point2D
    The end point of the ray segment.
length : float
    The length of the ray segment in meters.
travel_time : float
    The travel time along the ray segment in seconds.
)""";
    segment.def("__copy__", [](const Segment2D &self)
    {   
        return Segment2D(self);
    }); 
    segment.def("set_start_and_end_point",
                &Segment2D::setStartAndEndPoint,
                "Sets the starting and ending point of the ray segment.");
    segment.def("clear",
                &Segment2D::clear,
                "Resets the class.");
    segment.def_property_readonly("start_point",
                                  &Segment2D::getStartPoint);
    segment.def_property_readonly("end_point",
                                  &Segment2D::getEndPoint);
    segment.def_property_readonly("length",
                                  &Segment2D::getLength);
    segment.def_property_readonly("travel_time",
                                  &Segment2D::getTravelTime);
    segment.def_property("velocity",
                         &Segment2D::getVelocity,
                         &Segment2D::setVelocity);
    segment.def_property("velocity_model_cell_index",
                         &Segment2D::getVelocityModelCellIndex,
                         &Segment2D::setVelocityModelCellIndex);
}
