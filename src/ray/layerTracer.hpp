#ifndef PRIVATE_LAYER_TRACER_HPP
#define PRIVATE_LAYER_TRACER_HPP
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#ifndef NDEBUG
#include <cassert>
#endif
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"

namespace
{

enum class ReturnCode
{
    UnderShot =-1,
    Hit = 0,
    OverShot =+1,
    RayDoesNotTurn,
    RayTurnsTooEarly
};

struct Segment
{
    double slowness{0};
    double x0{0};
    double z0{0};
    double x1{0};
    double z1{0};
    int layer{0};
    void reverse()
    {
        std::swap(x0, x1);
        std::swap(z0, z1);
    }
    [[nodiscard]] EikonalXX::Ray::Segment2D toSegment() const
    {
        EikonalXX::Ray::Point2D startPoint{x0, z0};
        EikonalXX::Ray::Point2D endPoint{x1, z1};
        EikonalXX::Ray::Segment2D segment;
        segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
        segment.setSlowness(slowness);
        segment.setVelocityModelCellIndex(layer);
        return segment;
    }
    [[nodiscard]] double getTravelTime() const
    {
        return std::hypot(x1 - x0, z1 - z0)*slowness;
    }
};

/// Reverses a segment
void reverseSegments(std::vector<::Segment> *segments)
{
    std::reverse(segments->begin(), segments->end());
    for (auto &segment : *segments)
    {
        segment.reverse();
    }
}


[[nodiscard]] int getLayer(const double depth,
                           const std::vector<double> &interfaces,
                           const bool isAugmented = true)
{    
    if (interfaces.size() == 1){return 0;}
#ifndef NDEBUG
    if (isAugmented){assert(interfaces.size() > 1);}
#endif
    auto nLayers = static_cast<int> (interfaces.size());
    if (isAugmented){nLayers = nLayers - 1;}
    int layer = 0; 
    if (depth < interfaces.front())
    {
        layer = 0;
    }
    else if (depth >= interfaces.back())
    {
        layer = nLayers - 1;
    }
    else
    {
        layer = std::distance(
                   interfaces.begin(),
                   std::upper_bound(interfaces.begin(),
                                    interfaces.begin() + nLayers,
                                    depth)) - 1; 
#ifndef NDEBUG
        assert(depth >= interfaces.at(layer) &&
               depth < interfaces.at(layer + 1)); 
#endif
    }
#ifndef NDEBUG
    assert(layer >= 0);
    if (nLayers > 1){assert(layer < nLayers - 1);} 
#endif
    return layer;
}   

/// @brief Computes the critical angle from the slownesses.
/// @param[in] s0   The slowness in the current layer in s/m.
/// @param[in] s1   The slowness in the transmission layer in s/m.
/// @result The critical angle in radians.
/// @note Cerveny, Page 121
[[nodiscard]] double computeCriticalAngle(const double s0,
                                          const double s1)
{
    if (s0 <= s1){return M_PI_2;} // Use max transmission angle
    // arcsin(v0/v1) = arcsin( (1/s0)/(1/s1) ) = arcsin(s1/s0)
    return std::asin(s1/s0);
}

/// @brief Given the incidence angle and the slowness in this layer and the
///        next layer, this computes the corresponding transmission angle.
/// @param[in] incidenceAngle  The angle of incidence in radians.
/// @param[in] s0              The slowness in the current layer in s/m.
/// @param[in] s1              The slowness in the transmission laye rin s/m.
/// @result The transmission angle in radians.
/// @note Cerveny pg 42
[[nodiscard]] double computeTransmissionAngle(const double incidenceAngle,
                                              const double s0,
                                              const double s1)
{
#ifndef NDEBUG
    assert(incidenceAngle >= 0 && incidenceAngle <= M_PI);
#endif
    // sin(i0)/v0 = sin(i1)/v1 -> i1 = arcsin((s0/s1)*sin(i0))
    auto asinArgument = (s0/s1)*std::sin(incidenceAngle);
    if (asinArgument >= 0 && asinArgument <= 1)
    {
        return std::asin(asinArgument);
    }
    return -incidenceAngle; // Post critical -> total internal reflection
}

/// @result Defines the ray convergence
[[nodiscard]]
ReturnCode checkRayConvergence(const std::vector<::Segment> &segments,
                               const double stationOffset,
                               const double tolerance = 1)
{
    if (std::abs(segments.back().x1 - stationOffset) < tolerance)
    {
        return ReturnCode::Hit;
    }
    else
    {
        if (segments.back().x1 > stationOffset)
        {
            return ReturnCode::OverShot;
        }
        else
        {
            return ReturnCode::UnderShot;
        }
    }
#ifndef NDEBUG
    assert(false);
#endif
}

/// @brief Traces a ray strictly down through a stack of layers.
ReturnCode traceVerticalReflectionDown(const std::vector<double> &interfaces,
                                       const std::vector<double> &slownesses,
                                       const int sourceLayer,
                                       const int endLayer,
                                       const double sourceDepth,
                                       const int stationLayer,
                                       const double stationDepth,
                                       const double stationOffset,
                                       std::vector<::Segment> *segments,
                                       const double tolerance = 1)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1;
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(interfaces.size() == slownesses.size());
    assert(sourceLayer >= 0 && sourceLayer <= endLayer);
    assert(stationLayer >= 0);
    if (stationDepth > sourceDepth){assert(stationLayer <= endLayer);}
    assert(endLayer < nLayers - 1); // Can't bounce from whole-space
    assert(stationOffset >= 0); 
    assert(sourceDepth >= interfaces[sourceLayer] &&
           sourceDepth <  interfaces[sourceLayer + 1]);
    assert(stationDepth >= interfaces[stationLayer] &&
           stationDepth <  interfaces[stationLayer + 1]);
#endif
    // Same layer -> just get out of here
    constexpr double x0{0};
    constexpr double x1{0};
    double z0{sourceDepth};
    if (sourceLayer == stationLayer && sourceLayer == endLayer)
    {
        ::Segment segment{slownesses[sourceLayer],
                          x0, z0, x1, stationDepth,
                          sourceLayer};
        segments->push_back(segment);
        return ::checkRayConvergence(*segments, stationOffset, tolerance);
    } 
    if (segments->capacity() < 2*interfaces.size() + 1)
    {
        segments->reserve(2*interfaces.size() + 1);
    }
    // Trace from source down
    for (int i = sourceLayer; i <= endLayer; ++i)
    {
        double z1 = interfaces[i + 1]; 
        ::Segment segment{slownesses[i],
                          x0, z0, x1, z1, 
                          i};
        segments->push_back(segment);
        // Update
        z0 = z1;
    }
    // Unwind this thing back to the station 
    for (int i = endLayer; i >= stationLayer; --i)
    {
        double z0 = interfaces[i];
        if (i == stationLayer)
        {
            z0 = stationDepth;
        }
        ::Segment segment{slownesses[i],
                          x0,
                          interfaces[i + 1],
                          x1,
                          z0,
                          i};
        segments->push_back(segment);
    }
    // Check convergence
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// @brief Traces a symmetric ray path from a source, to a max
///        interface, then back up to a station at the same depth.
/// @param[in] interfaces    The (augmented) interfaces.  This indicates
///                          the depth of the top of each layer in meters
///                          which increases nadir.
/// @param[in] slownesses    The (agumented) slownesses in s/m in each
///                          layer.
/// @param[in] takeOffAngle  The take-off angle in degrees.
/// @param[in] sourceLayer   The layer containing the source.
/// @param[in] endLayer      The layer to which to trace.
/// @param[in] sourceDepth   The source depth in meters.
/// @param[in] stationOffset The source-receiver offset in meters.
/// @param[out] segments     The segments comprising the ray path.
/// @param[in] tolerance     A  hit is defined to be within this tolerance.
/// @result An indicator describing the ray convergence.
[[nodiscard]]
ReturnCode traceDown(const std::vector<double> &interfaces,
                     const std::vector<double> &slownesses,
                     const double takeOffAngle,
                     const int sourceLayer,
                     const int endLayer,
                     const double sourceDepth,
                     const double stationOffset,
                     std::vector<::Segment> *segments,
                     const double tolerance = 1)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1;
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(interfaces.size() == slownesses.size());
    assert(sourceLayer >= 0 && sourceLayer <= endLayer);
    assert(endLayer < nLayers);
    assert(takeOffAngle >= 0 && takeOffAngle < 90);
    assert(stationOffset >= 0);
    assert(sourceDepth >= interfaces[sourceLayer] &&
           sourceDepth <  interfaces[sourceLayer + 1]);
#endif
    if (segments->capacity() < 2*interfaces.size() + 1)
    {   
        segments->reserve(2*interfaces.size() + 1);
    } 
    double x0{0};
    double z0{sourceDepth};
    double currentAngle = takeOffAngle*(M_PI/180);
    for (int i = sourceLayer; i <= endLayer; ++i)
    {
        // Draw the ray
        double z1 = interfaces[i + 1];
        double x1 = x0 + (z1 - z0)*std::tan(currentAngle);
        ::Segment segment{slownesses[i],
                          x0, z0, x1, z1,
                          i};
        segments->push_back(segment);
        // Turning?
        auto criticalAngle 
            = computeCriticalAngle(slownesses[i],
                                   slownesses[i + 1]);
        if (currentAngle > criticalAngle)
        {
            if (i < endLayer - 1){return ReturnCode::RayTurnsTooEarly;}
            // Add the last segment by hand.  Compute the distance to the
            // half offset.
            auto dxHalf = 0.5*stationOffset - x1;
            // Now travel the mirror'ing critically refracted half
            // - i.e., 2*dxHalf
            ::Segment segment{slownesses[i + 1],
                              x1, z1, x1 + 2*dxHalf, z1,
                              i + 1};
            segments->push_back(segment);
            break; 
        }
        auto transmissionAngle
            = computeTransmissionAngle(currentAngle,
                                       slownesses[i],
                                       slownesses[i + 1]);
#ifndef NDEBUG
        assert(transmissionAngle >= 0);
#endif
        // Update
        x0 = x1;
        z0 = z1;
        currentAngle = transmissionAngle;
    }
    // Unwind this thing   
    x0 = segments->back().x1;
    for (int i = endLayer - sourceLayer; i >= 0; --i)
    {
        auto dx = segments->at(i).x1 - segments->at(i).x0;
        ::Segment segment{segments->at(i).slowness,
                          x0,
                          segments->at(i).z1,
                          x0 + dx,
                          segments->at(i).z0,
                          segments->at(i).layer};
        segments->push_back(segment);
        x0 = x0 + dx;
    }
    // Check convergence
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// Traces from a direct ray upwards
[[nodiscard]]
ReturnCode traceDirect(const std::vector<double> &interfaces,
                       const std::vector<double> &slownesses,
                       const double takeOffAngle,
                       const int sourceLayer,
                       const double sourceDepth,
                       const int stationLayer,
                       const double stationOffset,
                       const double stationDepth,
                       std::vector<::Segment> *segments,
                       const double tolerance = 1)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1; 
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(stationLayer <= sourceLayer);
    assert(takeOffAngle > 90 && takeOffAngle <= 180);
    assert(sourceLayer >= 0 && sourceLayer < nLayers);
    assert(stationOffset >= 0);
    assert(sourceDepth < stationDepth);
    assert(sourceDepth >= interfaces[sourceLayer] &&
           sourceDepth <  interfaces[sourceLayer + 1]);
    assert(stationDepth >= interfaces[stationLayer] &&
           stationDepth <  interfaces[stationLayer + 1]);
#endif
    segments->reserve(sourceLayer - stationLayer + 1);
    // Simplify geometry
    double currentAngle = std::abs(180 - takeOffAngle)*(M_PI/180);
    if (takeOffAngle == 180){currentAngle = 0;}
    double x0{0};
    double z1{sourceDepth};
    for (int layer = sourceLayer; layer >= stationLayer; --layer) 
    {
        double z0 = interfaces[layer];
        if (layer == stationLayer){z0 = stationDepth;}
        double x1 = x0 + (z1 - z0)*std::tan(currentAngle);
        if (takeOffAngle == 180){x1 = x0;}
        ::Segment segment{slownesses[layer],
                          x0, z0, x1, z1,
                          layer};
        segments->push_back(segment);
        if (layer == stationLayer){break;}
        double transmissionAngle
            = computeTransmissionAngle(currentAngle,
                                       slownesses[layer],
                                       slownesses.at(layer - 1));
        //double criticalAngle = computeCriticalAngle(slownesses[layer],
        //                                            slownesses[layer - 1]);
        // Update
        x0 = x1;
        z1 = z0;
        currentAngle = transmissionAngle;
    }
    // Check convergence
    return checkRayConvergence(*segments, stationOffset, tolerance);
}

/// 
[[nodiscard]]
ReturnCode traceDownThenUp(const std::vector<double> &interfaces,
                           const std::vector<double> &slownesses,
                           const double takeOffAngleIn,
                           const int sourceLayerIn,
                           const double sourceDepthIn,
                           const int stationLayerIn,
                           const double stationOffset,
                           const double stationDepthIn,
                           const int endLayer,
                           std::vector<::Segment> *segments,
                           const double tolerance = 1)
{
    auto takeOffAngle = takeOffAngleIn;
#ifndef NDEBUG
    assert(takeOffAngle >= 0 && takeOffAngle < 90);
    assert(segments != nullptr);
#endif
    segments->clear();
    if (std::abs(stationDepthIn - sourceDepthIn) < 1.e-10)
    {
        return ::traceDown(interfaces,
                           slownesses,
                           takeOffAngle,
                           sourceLayerIn,
                           endLayer,
                           sourceDepthIn,
                           stationOffset,
                           segments,
                           tolerance);
    }
    // Reverse segments? 
    bool reverseSegments{false};
    if (stationDepthIn < sourceDepthIn){reverseSegments = true;}
    // General ray tracing
    auto sourceLayer = sourceLayerIn;
    auto sourceDepth = sourceDepthIn;
    auto stationLayer = stationLayerIn;
    auto stationDepth = stationDepthIn;
    // Trace out the top segments
    auto takeOffAngleUp = takeOffAngle + 90;
    std::vector<::Segment> upSegments;
    auto result = ::traceDirect(interfaces,
                                slownesses,
                                takeOffAngleUp,
                                sourceLayer,
                                sourceDepth,
                                stationLayer,
                                stationOffset,
                                stationDepth,
                                &upSegments,
                                tolerance);
    if (result != ReturnCode::Hit &&
        result != ReturnCode::UnderShot &&
        result != ReturnCode::OverShot)
    {
        return result;
    }
    // Trace the bottom segments from source to `source'
    auto offsetUp = upSegments.back().x1 - upSegments.front().x0;
    std::vector<::Segment> downSegments;
    result = ::traceDownThenUp(interfaces,
                               slownesses,
                               takeOffAngle,
                               sourceLayer,
                               sourceDepth,
                               sourceLayer, 
                               offsetUp,
                               sourceDepth,
                               endLayer,
                               &downSegments,
                               tolerance);
    if (result != ReturnCode::Hit &&
        result != ReturnCode::UnderShot &&
        result != ReturnCode::OverShot)
    {
        return result;
    }
    auto downOffset = downSegments.back().x1;
    // Put the ray path together 
    segments->resize(downSegments.size() + upSegments.size());
    std::copy(downSegments.begin(), downSegments.end(),
              segments->begin());
    int j = static_cast<int> (downSegments.size());
    for (auto &upSegment : upSegments)
    {
        upSegment.x0 = upSegment.x0 + downOffset;
        upSegment.x1 = upSegment.x1 + downOffset;
        segments->at(j) = std::move(upSegment);
    }
    if (reverseSegments){::reverseSegments(segments);}
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// @brief Converts the 1D velocity stack to a slowness stack. 
[[nodiscard]]
std::vector<double> toSlownessVector(const std::vector<double> &velocities)
{
    std::vector<double> slownesses(velocities.size());
    std::transform(velocities.begin(), velocities.end(), slownesses.begin(),
                   [&](const double velocity)
                   {
#ifndef NDEBUG
                       assert(velocity >= 0);
#endif
                       return 1./velocity;
                   });
    return slownesses;
}

/// @brief Modifies the velocity vector for sfae use with the above routines.
[[nodiscard]] 
std::vector<double> augmentVelocityVector(const std::vector<double> &x)
{
    auto y = x;
    if (!y.empty()){y.push_back(x.back());}
    return y;
}

/// @brief Modifies the interfaces vector for safe use with the above routines.
[[nodiscard]]
std::vector<double> augmentInterfacesVector(const std::vector<double> &x)
{
    auto y = x;
    y.push_back(std::numeric_limits<double>::max());
    return y;
}

/// @brief Converts the local ray path to the library variant.
[[nodiscard]]
EikonalXX::Ray::Path2D
    toRayPath(const std::vector<::Segment> &segments)
{
    if (segments.empty()){throw std::runtime_error("Segments is empty");}
    std::vector<EikonalXX::Ray::Segment2D> raySegments;
    raySegments.reserve(segments.size());
    for (const auto &segment : segments)
    {
        raySegments.push_back(segment.toSegment());
    }
    EikonalXX::Ray::Path2D rayPath;
    rayPath.set(std::move(raySegments));
    return rayPath;
}

/// Traces a direct wave when the source and receiver are in the smae layer.
[[nodiscard]]
ReturnCode
    traceDirectSameLayer(const std::vector<double> &interfaces,
                         const std::vector<double> &slownesses,
                         const int sourceLayer,
                         const double sourceDepth,
                         const int stationLayer,
                         const double stationOffset,
                         const double stationDepth,
                         EikonalXX::Ray::Path2D *path)
{
    path->clear();
#ifndef NDEBUG
    assert(sourceLayer == stationLayer);
    assert(!slownesses.empty());
    assert(slownesses.size() == interfaces.size());
    assert(path != nullptr);
#endif
    EikonalXX::Ray::Point2D startPoint{0, sourceDepth};
    EikonalXX::Ray::Point2D endPoint{stationOffset, stationDepth};
    EikonalXX::Ray::Segment2D segment;
    segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
    segment.setSlowness(slownesses.at(sourceLayer));
    segment.setVelocityModelCellIndex(sourceLayer);
    path->open();
    path->append(segment);
    path->close();
    return ReturnCode::Hit;
}

/// @brief Traces a direct ray in a whole space.
[[nodiscard]]
ReturnCode
    traceWholeSpace(const double slowness,
                    const double sourceDepth,
                    const double stationOffset,
                    const double stationDepth,
                    EikonalXX::Ray::Path2D *path)
{
#ifndef NDEBUG
    assert(slowness >= 0);
    assert(path != nullptr);
#endif
    const std::vector<double> interfaces{std::numeric_limits<double>::min(),
                                         std::numeric_limits<double>::max()};
    const std::vector<double> slownesses{slowness, slowness};
    return ::traceDirectSameLayer(interfaces, slownesses,
                                  0, sourceDepth,
                                  0, stationOffset, stationDepth,
                                  path);
}

}
 
/*
int main()
{
    const std::vector<double> interfaces{-4500,   50, 15600, 26500, 40500};
    const std::vector<double> velocities{3500, 5900,  6400,  7500,  7900};
    const std::vector<double> offsets{15000, 20000, 30000, 40000, 50000, 60000, 80000, 100000, 120000};
    const std::vector<double> eikonalTimes{3.4752889922,
                                           4.32274661932,
                                           6.01766187355,
                                           7.71257712779,
                                           9.40749238203,
                                           11.1024076363,
                                           14.4922381447,
                                           17.8820686532,
                                           21.2718991617};
    auto augmentedInterfaces = ::augmentInterfacesVector(interfaces);
    auto augmentedSlownesses
        = ::toSlownessVector(::augmentVelocityVector(velocities));
    double takeOffAngle{39.2777};
    double sourceDepth{-2000};
    auto stationDepth = sourceDepth;
    auto sourceLayer = ::getLayer(sourceDepth, augmentedInterfaces, true);
    auto stationLayer = ::getLayer(stationDepth, augmentedInterfaces, true);
    for (size_t iOffset = 0; iOffset < offsets.size(); ++iOffset)
    {
        std::vector<::Segment> segments;
        auto result = ::traceDownThenUp(augmentedInterfaces,
                                        augmentedSlownesses,
                                        takeOffAngle,
                                        sourceLayer,
                                        sourceDepth,
                                        stationLayer,
                                        offsets[iOffset],
                                        stationDepth,
                                        1,
                                        &segments,
                                        1.0);
        if (result == ReturnCode::Hit || result == ReturnCode::UnderShot || result == ReturnCode::OverShot)
        {
            double ttTime = 0;
            for (auto &segment : segments)
            {
                ttTime = ttTime + segment.getTravelTime();
                //std::cout << segment.x0 << "," << segment.z0 << ","
                //          << segment.x1 << "," << segment.z1 << ","
                //          << 1./segment.slowness << std::endl;
            }
            //std::cout << std::setprecision(12) << offsets[iOffset] << " " << ttTime << " " << eikonalTimes[iOffset] << " " << ttTime - eikonalTimes[iOffset] << std::endl;
        }
    }
    // Trace some direct waves
    std::vector<::Segment> upSegments;
    sourceDepth = 30000;
    sourceLayer = ::getLayer(sourceDepth, augmentedInterfaces, true);
    auto result = ::traceDirect(augmentedInterfaces,
                              augmentedSlownesses,
                              150,
                              sourceLayer,
                              sourceDepth,
                              stationLayer,
                              12410, //offsets[0],
                              stationDepth,
                              &upSegments,
                              1.01);
    double ttime = 0;
    for (const auto &segment : upSegments)
    {
        //std::cout << "segment time: " << segment.getTravelTime() << std::endl;
        ttime = ttime + segment.getTravelTime();
    }
    std::cout << ttime << std::endl;
    return EXIT_SUCCESS;
}
*/

#endif