#ifndef PYEIKONALXX_RAY_2D_HPP
#define PYEIKONALXX_RAY_2D_HPP
#include <memory>
#include <vector>
#include <pybind11/pybind11.h>
namespace EikonalXX::Ray
{
class Point2D;
class Segment2D;
}
namespace PEikonalXX::Ray
{
class Point2D
{
public:
    Point2D();
    Point2D(const Point2D &point);
    Point2D(const EikonalXX::Ray::Point2D &point);
    Point2D(Point2D &&point) noexcept;
    Point2D& operator=(const Point2D &station);
    Point2D& operator=(const EikonalXX::Ray::Point2D &point);
    Point2D& operator=(Point2D &&point) noexcept;
    void setPositionInX(double x) noexcept;
    [[nodiscard]] double getPositionInX() const;
    [[nodiscard]] bool havePositionInX() const noexcept;
    void setPositionInZ(double z) noexcept;
    [[nodiscard]] double getPositionInZ() const;
    [[nodiscard]] bool havePositionInZ() const noexcept;
    [[nodiscard]] EikonalXX::Ray::Point2D& getNativeClassReference() const;
    ~Point2D();
    void clear() noexcept;
private:
    std::unique_ptr<EikonalXX::Ray::Point2D> pImpl;
};
class Segment2D
{
public:
    Segment2D();
    Segment2D(const Segment2D &segment);
    Segment2D(const EikonalXX::Ray::Segment2D &segment);
    Segment2D(Segment2D &&segment) noexcept;
    Segment2D& operator=(const Segment2D &segment);
    Segment2D& operator=(const EikonalXX::Ray::Segment2D &segment);
    Segment2D& operator=(Segment2D &&segment) noexcept;
    void setStartAndEndPoint(const std::pair<Point2D, Point2D> &startAndEndPoint);
    [[nodiscard]] Point2D getStartPoint() const;
    [[nodiscard]] Point2D getEndPoint() const;
    [[nodiscard]] bool haveStartAndEndPoint() const noexcept;
    [[nodiscard]] double getLength() const;
    void setVelocity(double velocity); 
    void setSlowness(double slowness);
    [[nodiscard]] double getVelocity() const;
    [[nodiscard]] double getSlowness() const;
    [[nodiscard]] bool haveVelocity() const noexcept;
    [[nodiscard]] double getTravelTime() const;
    void setVelocityModelCellIndex(int index);
    [[nodiscard]] int getVelocityModelCellIndex() const;
    [[nodiscard]] bool haveVelocityModelCellIndex() const noexcept;
    [[nodiscard]] EikonalXX::Ray::Segment2D& getNativeClassReference() const;
    ~Segment2D();
    void clear() noexcept;
private:
    std::unique_ptr<EikonalXX::Ray::Segment2D> pImpl;
};
void initializeRay2D(pybind11::module &m);
}
#endif
