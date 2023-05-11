#include <cmath>
#include <vector>
#include <limits>
#include "eikonalxx/ray/path3d.hpp"
#include "eikonalxx/ray/segment3d.hpp"
#include "eikonalxx/ray/point3d.hpp"
#include "length.hpp"

using namespace EikonalXX::Ray;

namespace
{
void checkSegment(const Segment3D &segment)
{
    if (!segment.haveStartAndEndPoint())
    {
        throw std::invalid_argument("Segment start/end point not set");
    }
    if (!segment.haveVelocity())
    {
        throw std::invalid_argument("Segment velocity not set");
    }
}
}

class Path3D::Path3DImpl
{
public:
    void integrate()
    {
        double length = 0;
        double travelTime = 0;
        for (const auto &segment : mSegments)
        {
            travelTime = travelTime + segment.getTravelTime();
            length = length + segment.getLength();
        }
        mTravelTime = travelTime;
        mLength = length;
    }
    std::vector<Segment3D> mSegments;
    std::vector<Segment3D> mSegmentsConstruction;
    double mTravelTime{0};
    double mLength{0};
    bool mOpened{false};
};

/// Constructor
Path3D::Path3D() :
    pImpl(std::make_unique<Path3DImpl> ())
{
}

/// Copy constructor
Path3D::Path3D(const Path3D &path)
{
    *this = path;
}

/// Move constructor
Path3D::Path3D(Path3D &&path) noexcept
{
    *this = std::move(path);
}

/// Reset class
void Path3D::clear() noexcept
{
    pImpl = std::make_unique<Path3DImpl> ();
}

/// Destructor
Path3D::~Path3D() = default;

/// Copy assignment
Path3D& Path3D::operator=(const Path3D &path)
{
    if (&path == this){return *this;}
    pImpl = std::make_unique<Path3DImpl> (*path.pImpl);
    return *this;
}

/// Move constructor
Path3D& Path3D::operator=(Path3D &&path) noexcept
{
    if (&path == this){return *this;}
    pImpl = std::move(path.pImpl);
    return *this;
}

/// Open the ray path
void Path3D::open()
{
    pImpl->mSegmentsConstruction.clear();
    pImpl->mSegmentsConstruction.reserve(
        std::max(128, static_cast<int> (pImpl->mSegments.size())));
    pImpl->mOpened = true;
}

bool Path3D::isOpen() const noexcept
{
    return pImpl->mOpened;
}

/// Append to the ray path
void Path3D::append(const Segment3D &segment)
{
    auto temporarySegment = segment;
    append(std::move(temporarySegment));
}

void Path3D::append(Segment3D &&segment)
{
    ::checkSegment(segment);
    if (pImpl->mSegmentsConstruction.empty())
    {
        pImpl->mSegmentsConstruction.push_back(std::move(segment));
    }
    else
    {
        auto point0 = pImpl->mSegmentsConstruction.back().getEndPoint();
        auto point1 = segment.getStartPoint();
        auto distance = ::computeLength(point0, point1);
        if (distance > std::numeric_limits<double>::epsilon()*100)
        {
            throw std::invalid_argument(
               "Segment does no start at last segments's end point"); 
        }
        pImpl->mSegmentsConstruction.push_back(std::move(segment));
    }
}

void Path3D::set(const std::vector<Segment3D> &segments)
{
    auto temporarySegments = segments;
    set(std::move(temporarySegments));
}

void Path3D::set(std::vector<Segment3D> &&segments)
{
    for (int i = 0; i < static_cast<int> (segments.size()); ++i)
    {
        ::checkSegment(segments[i]);
        if (i > 0)
        {
            auto distance = ::computeLength(segments[i-1].getEndPoint(),
                                            segments[i].getStartPoint());
            if (distance > std::numeric_limits<double>::epsilon()*100)
            {
                throw std::invalid_argument("Segment " + std::to_string(i)
                                          + " does no start at last segments "
                                          + std::to_string(i - 1)
                                          + " end point");
            }
        }
    }
    pImpl->mSegments = std::move(segments);
    pImpl->mSegmentsConstruction.clear();
    pImpl->mOpened = false;             
    // sum along the path
    pImpl->integrate();
}

/// Close the ray path
void Path3D::close()
{
    if (!isOpen()){throw std::runtime_error("Path not open");}
    // Do some tap dancing and release memory
    pImpl->mSegments = std::move(pImpl->mSegmentsConstruction);
    pImpl->mSegmentsConstruction = std::vector<Segment3D> ();
    pImpl->mOpened = false;
    // Sum along the path
    pImpl->integrate();
}

size_t Path3D::size() const noexcept
{
    return static_cast<int> (pImpl->mSegments.size());
}

/// Travel time
double Path3D::getTravelTime() const
{
    return pImpl->mTravelTime;
}

/// Length
double Path3D::getLength() const
{
    return pImpl->mLength;
}

/// Iterators
Path3D::iterator Path3D::begin()
{
    return pImpl->mSegments.begin();
}

constexpr Path3D::const_iterator Path3D::begin() const noexcept
{
    return pImpl->mSegments.begin();
}

constexpr Path3D::const_iterator Path3D::cbegin() const noexcept
{
    return pImpl->mSegments.cbegin();
}

Path3D::iterator Path3D::end()
{
    return pImpl->mSegments.end();
}

constexpr Path3D::const_iterator Path3D::end() const noexcept
{
    return pImpl->mSegments.end();
}

constexpr Path3D::const_iterator Path3D::cend() const noexcept
{
    return pImpl->mSegments.cend();
}

Segment3D& Path3D::at(const size_t index)
{
    return pImpl->mSegments.at(index);
}

const Segment3D& Path3D::at(const size_t index) const
{
    return pImpl->mSegments.at(index);
}

