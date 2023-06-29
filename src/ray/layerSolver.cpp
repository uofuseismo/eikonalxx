#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include "eikonalxx/ray/layerSolver.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/point2d.hpp"

using namespace EikonalXX::Ray;

namespace
{

enum class RayType
{
    Direct,
    Reflected,
    Refracted
};

enum class ReturnCode
{
    Undershot,    /*!< The ray exited the of the medium before the station offset. */
    Overshot,     /*!< The ray exited the medium after the station offset. */
    ExitedBottom, /*!< The ray exited the bottom of the medium. */
    ExitedRight,  /*!< The ray exited the right side of the medium. */
    Hit,          /*!< The ray hit its intended target to within some tolerance. */
    RayDoesNotTurn,
    RayTurnsTooEarly
};

struct Segment
{
    [[nodiscard]] Segment2D toSegment() const
    {
        Segment2D result;
        Point2D point0{x0, z0};
        Point2D point1{x1, z1};
        result.setStartAndEndPoint(std::pair {point0, point1});
        result.setSlowness(slowness);
        result.setVelocityModelCellIndex(layerIndex);
        return result;
    } 
    double x0;
    double z0;
    double x1;
    double z1;
    double slowness;
    int layerIndex;
};

/// Traces down
/// ------------- Interface 0
/// x          x    Slowness 0       
/// ------------- Interface 1
///   \      /      Slowness 1
///-------------- Interface 2
///     \__/        Slowness 2
///-------------- Interace 3
///                 Slowness 3
ReturnCode levelTraceDown(const double depth,
                          const double offset,
                          const double takeOffAngle,
                          const int startIndex,
                          const int endIndex,
                          const std::vector<double> &interfaces,
                          const std::vector<double> &slownesses,
                          Path2D *path,
                          const double tolerance = 10,
                          const bool returnToSurface = true)
{
    path->clear();
    // Just deal with same layer problem
    if (startIndex == endIndex)
    {
#ifndef NDEBUG
        std::cerr << __func__ << " startIndex == endIndex!" << std::endl;
#endif
        Point2D point0{0, depth};
        Point2D point1{offset, depth}; 
        Segment2D segment;
        segment.setStartAndEndPoint(std::pair {point0, point1});
        segment.setSlowness(slownesses.at(startIndex));
        segment.setVelocityModelCellIndex(startIndex);
        path->clear();
        path->open();
        path->append(segment);
        path->close();
        return ReturnCode::Hit;
    }
    // This will be problematic
    if (takeOffAngle >= 90)
    {
        std::cerr << "Ray take-off angle must be less than 90" << std::endl;
        return ReturnCode::RayDoesNotTurn;
    }
    auto nLayers = static_cast<int> (interfaces.size());
    constexpr double degreesToRadians{M_PI/180};
    auto takeOffAngleRadians = takeOffAngle*degreesToRadians;
    //double p = std::sin(takeOffAngleRadians)*slownesses.at(startIndex);
    auto halfOffset = 0.5*offset;
    double x0{0};
    double z0{depth};
    double z1 = interfaces[startIndex + 1];
    double dz = z1 - depth;
#ifndef NDEBUG
    assert(dz >= 0);
#endif
    double x1 = 0;
    if (dz > 0)
    {
        x1 = dz*std::tan(takeOffAngleRadians);///std::cos(takeOffAngleRadians);
    }
    // Trace out first layer
    std::vector<::Segment> segments;
    segments.reserve(2*nLayers);
    ::Segment firstSegment{x0, z0, x1, z1,
                           slownesses.at(startIndex), startIndex};
    segments.push_back(std::move(firstSegment));
    // Now trace through the stack 
    bool criticallyRefracted = false;
    auto currentAngle = takeOffAngleRadians;
    for (int layer = startIndex + 1; layer <= endIndex; ++layer)
    {
        auto s0 = slownesses.at(layer - 1);
        auto s1 = slownesses.at(layer);
        auto s0s1 = s0/s1; // V_1/V_0 = (1/s_1)/(1/s_0) = s_0/s_1
        // Update positions
        x0 = x1;
        z0 = z1;
        z1 = interfaces.at(layer + 1);
        dz = z1 - z0;
        x1 = dz*std::tan(currentAngle);///std::cos(currentAngle);
        // Swing and a miss 
        if (x1 > halfOffset)
        {
            return ReturnCode::RayDoesNotTurn;
        } 
        ::Segment segment{x0, z0, x1, z1, 
                          s1, layer};
        segments.push_back(std::move(segment));
        // Note the critical angle
        double criticalAngle{M_PI_2}; 
        if (s0s1 > 1){criticalAngle = std::asin(s1/s0);} 
        // Critical refraction.  Simply move the segment
        // along the interface to the half-way point.
        if (currentAngle > criticalAngle)
        {
            if (layer != endIndex)
            {
                //std::cout << "Ray critically refracted before end" << std::endl;
                return ReturnCode::RayTurnsTooEarly;
            }
            // Went too far.  Just bounce back to the surface.
            if (x1 >= halfOffset){break;}
            //std::cout << "Ray critically refracted" << std::endl;
            auto dOffset2 = 0.5*offset - x1;
            ::Segment segment{x1, z1, 
                              x1 + 2*dOffset2, z1,
                              slownesses.at(layer + 1), layer + 1};
            segments.push_back(segment);
            criticallyRefracted = true;
            break;
        }
        // Update angle with Snell's law
        auto asinArgument = std::sin(currentAngle)*s0s1;
#ifndef NDEBUG
        assert(asinArgument >=-1 && asinArgument <= 1);
#endif
        //std::cout << "layer-1, layer, current angle, critical angle: " << layer - 1 << " " << layer << " " << currentAngle/degreesToRadians << "," << criticalAngle/degreesToRadians << std::endl;
        currentAngle = std::asin(asinArgument);
    }
    // Take the ray back to the `surface'
    if (returnToSurface)
    {
        auto iStart = static_cast<int> (segments.size()) - 1;
        if (criticallyRefracted)
        {
            iStart = iStart - 1;
        }
        x0 = segments.back().x1;
        for (int i = iStart; i >= 0; --i)
        {
            auto dx = segments[i].x1 - segments[i].x0;
#ifndef NDEBUG
            assert(dx >= 0);
#endif
            x1 = x0 + dx;
            ::Segment newSegment{x0, segments[i].z1,
                                 x1, segments[i].z0,
                                 segments[i].slowness,
                                 segments[i].layerIndex};
            segments.push_back(std::move(newSegment));
            // Update x
            x0 = x1;
        }
    }
    // Build the ray path
    path->clear();
    path->open();
    for (auto &segment : segments)
    {
        path->append(segment.toSegment());
    }
    path->close();
    if (std::abs(x1 - offset) < tolerance)
    {
        return ReturnCode::Hit;
    }
    else
    {
        if (x1 < offset)
        {
            return ReturnCode::Undershot;
        }
    }
    return ReturnCode::Overshot;
}

ReturnCode traceDown(const int stationLayer,
                     const int sourceLayer,
                     const double stationDepth,
                     const double sourceDepth,
                     const double offset,
                     const double takeOffAngle,
                     const int startIndex,
                     const int endIndex,
                     const std::vector<double> &interfaces,
                     const std::vector<double> &slownesses,
                     Path2D *path,
                     const double tolerance = 10)
{
    // Ray parameter is constant so we can get the station incidence angle
    // and make that our take-off angle to trace down to the source depth.
    auto takeOffAngleRadians = takeOffAngle*(M_PI/180);
    double p = std::sin(takeOffAngleRadians)*slownesses.at(startIndex);
    auto stationSlowness = slownesses.at(stationLayer);
    auto asinArgument = p/stationSlowness;
    double incidenceAngle = M_PI_2;
    if (asinArgument >= -1 && asinArgument <= 1)
    {
        incidenceAngle = std::asin(p/stationSlowness);
    }
    else
    {
        std::cerr << "no clue" << std::endl;
    }        
    // Trace from station to source 
    auto nTopLayers = sourceLayer - stationLayer;
    std::vector<double> temporaryInterfaces(nTopLayers);
    std::copy(interfaces.begin(), interfaces.begin() + nTopLayers,
              temporaryInterfaces.begin());
    temporaryInterfaces.back() = sourceDepth;
    std::vector<double> temporarySlownesses(nTopLayers);
    std::copy(slownesses.begin(), slownesses.begin() + nTopLayers,
              temporarySlownesses.begin());
    temporarySlownesses.back() = slownesses[sourceLayer];

    Path2D directPath;
    auto returnCode = levelTraceDown(stationDepth,
                                     0.5*offset, // If we go past this then we don't stand a chance
                                     incidenceAngle, // Now the takeoff angle
                                     stationLayer,
                                     temporaryInterfaces.size() - 1,
                                     temporaryInterfaces,
                                     temporarySlownesses,
                                     &directPath,
                                     tolerance,
                                     false);
}

void computeTakeOffAngles(
    const int n, const double theta0, const double theta1,
    std::vector<double> *takeOffAngles)
{
#ifndef NDEBUG
    assert(n > 0);
#endif
    if (static_cast<int> (takeOffAngles->size()) != n)
    {
        takeOffAngles->resize(n);
    }
    if (n == 1)
    {
        takeOffAngles->at(0) = 0.5*(theta0 + theta1);
    }
    auto dTheta = (theta1 - theta0)/(n - 1);
    auto *__restrict__ takeOffAnglePointer = takeOffAngles->data();
    for (int i = 0; i < n; ++i)
    {
        takeOffAnglePointer[i] = theta0 + dTheta*i;
    }
}

}

class LayerSolver::LayerSolverImpl
{
public:
    /// For when the source and receiver are in the same layer this computes
    /// the direct ray. 
    [[nodiscard]] Path2D computeDirectRayInSameLayer() const
    {
        Path2D rayPath;
#ifndef NDEBUG
        assert(mSourceLayer == mStationLayer);
#endif
        Segment2D segment;
        Point2D startPoint{0, mSourceDepth};
        Point2D endPoint{mStationOffset, mStationDepth};
        segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
        segment.setSlowness(mSlownessModel[mSourceLayer]);
        segment.setVelocityModelCellIndex(mSourceLayer);
        rayPath.open();
        rayPath.append(std::move(segment));
        rayPath.close();
        return rayPath;
    }
    [[nodiscard]] int getLayer(const double depth) const
    {
        auto nLayers = static_cast<int> (mInterfaces.size());
        int layer = 0;
        if (depth < mInterfaces.front())
        {
            layer = 0;
        }
        else if (depth >= mInterfaces.back())
        {
            layer = nLayers - 1;
        }
        else
        {
            layer = std::distance(
                       mInterfaces.begin(),
                       std::upper_bound(mInterfaces.begin(),
                                        mInterfaces.end(),
                                        depth)) - 1;
             //std::cout << "layer: " << mInterfaces[layer] << "," << depth << "," << mInterfaces[layer+1]<< std::endl;
#ifndef NDEBUG
            assert(depth >= mInterfaces.at(layer) &&
                   depth < mInterfaces.at(layer + 1));
#endif
        }
#ifndef NDEBUG
        assert(layer >= 0);
        if (nLayers > 1){assert(layer < nLayers - 1);}
#endif
        return layer;
    }
    /// @param[in] p   The ray parameter (sin(i)/v) where i is the take-off
    ///                angle measured positive up from nadir at the source
    ///                and v the velocity at the source.
    /// @note This cannot handle velocity inversions. 
    ReturnCode shoot(const double takeOffAngle)
    {
        Path2D path;
auto offset = mStationOffset;
        auto nLayers = static_cast<int> (mInterfaces.size());
        for (int endLayer = mSourceLayer + 1; endLayer < nLayers - 1; ++endLayer)
        {
            auto returnCode = ::levelTraceDown(mSourceDepth,
                                               offset,
                                               takeOffAngle,
                                               mSourceLayer,
                                               endLayer,
                                               mInterfaces,
                                               mSlownessModel,
                                               &path);
            if (returnCode == ReturnCode::Hit ||
                returnCode == ReturnCode::Undershot ||
                returnCode == ReturnCode::Overshot)
            {
                std::cout << std::endl;
                for (const auto &segment : path)
                {
                    std::cout << std::setprecision(12)
                              << segment.getStartPoint().getPositionInX() << ","
                              << segment.getStartPoint().getPositionInZ() << ","
                              << segment.getEndPoint().getPositionInX() << ","
                              << segment.getEndPoint().getPositionInZ() << ","
                              << 1./segment.getSlowness() << "," 
                              << 1./mSlownessModel[segment.getVelocityModelCellIndex()] << std::endl;
                }
                if (returnCode == ReturnCode::Hit)
                {
                    std::cout << "Takeoff angle, travel time: " << takeOffAngle << "," << path.getTravelTime() << std::endl;
                }
            }
        }
/*
getchar();

        if (takeOffAngle >= 180)
        {
            return ReturnCode::ExitedBottom;
        }
        constexpr int maxSegments{std::numeric_limits<int>::max()};
        constexpr double degreesToRadians{M_PI/180};
        constexpr double radiansToDegrees{180/M_PI};
        auto sourceSlowness = mSlownessModel[mSourceLayer]; 
        const double rayParameter
        {
            sourceSlowness*std::sin(takeOffAngle*degreesToRadians)
        };
        double x0 = 0;
        double z0 = mSourceDepth;
        auto nLayers = static_cast<int> (mSlownessModel.size());
        auto currentAngle = takeOffAngle;
        auto currentLayer = mSourceLayer;
        for (int kSegment = 0; kSegment < maxSegments; ++kSegment)
        {
            int nextLayerIncrement =+1; // Going down
            if (currentAngle > 90)
            {
                nextLayerIncrement = 0;
            }
            // Go to next interface
            double dz = 0;
            if (nextLayerIncrement == 1)
            {
                dz = mInterfaces[currentLayer + nextLayerIncrement] - z0;
            } 
            else
            {
                dz = z0 - mInterfaces[currentLayer + nextLayerIncrement];
            }
            double z1 = mInterfaces[currentLayer + nextLayerIncrement];
            // Laterally moving ray
            if (std::abs(currentAngle - 90) < 1.e-8)
            {
                Point2D endPoint{mMaximumOffset, z0}; 
                return ReturnCode::ExitedRight;
            } 
           
            double ds = dz/std::cos(currentAngle*degreesToRadians);
std::cout << std::setprecision(12) << dz << "," << ds << std::endl;
#ifndef NDEBUG
            assert(ds >= z1);
#endif 
            double x1 = x0 + std::sqrt(ds*ds - z1*z1);
            // Update the angle in the next layer with Snell's law 
            double slownessAbove = 0;
            double slownessBelow = 0;
            if (nextLayerIncrement == 1)
            {
                if (currentLayer == nLayers - 1)
                {
                    return ReturnCode::ExitedBottom;
                }
                slownessAbove = mSlownessModel[currentLayer];
                slownessBelow = mSlownessModel[currentLayer + 1];
            }
            else
            {
                if (currentLayer == 0)
                {
                    std::cout << "exited top" << std::endl;
                    getchar();
//                    return ReturnCode::ExitedTop;
                }
                slownessAbove = mSlownessModel[currentLayer - 1];
                slownessBelow = mSlownessModel[currentLayer];
            }
            // Next layer is faster than this layer - we're alright
            if (slownessAbove >= slownessBelow)
            {
                currentAngle = std::asin(rayParameter/slownessBelow)*radiansToDegrees;
            }
 std::cout << currentAngle << std::endl;
            // Ray exited bottom of medium
            // Update
            
getchar();
        }
*/
    }
    std::vector<Path2D> mRayPaths;
    std::vector<double> mSlownessModel;
    std::vector<double> mInterfaces;
    double mMaximumOffset{500000};
    double mStationDepth{0};
    double mSourceDepth{0};
    double mStationOffset{0};
    int mSourceLayer{0};
    int mStationLayer{0};
    bool mHaveStationOffsetAndDepth{false};
    bool mHaveSourceDepth{false};
    bool mHaveRayPaths{false};
};

/// Constructor
LayerSolver::LayerSolver() :
    pImpl(std::make_unique<LayerSolverImpl> ())
{
}

/*
/// Move constructor
LayerSolver::LayerSolve(LayerSolver &&solver) noexcept
{
    *this = std::move(solver);
}

/// Move assignment
LayerSolver& LayerSolver::operator=(LayerSolver &&solver) noexcept
{
    if (&solver == this){return *this;}
    pImpl = std::move(solver.pImpl);
    return *this;
}
*/

/// Reset class
void LayerSolver::clear() noexcept
{
    pImpl = std::make_unique<LayerSolverImpl> ();
}

/// Destructor
LayerSolver::~LayerSolver() = default;

/// Set the velocity model
void LayerSolver::setVelocityModel(const std::vector<double> &interfaces,
                                   const std::vector<double> &velocityModel)
{
    if (velocityModel.empty())
    {
        throw std::invalid_argument("Velocity model is empty");
    }
    if (interfaces.size() != velocityModel.size())
    {
        throw std::invalid_argument(
            "Velocity model size must equal interfaces size"
        );
    }
    for (const auto &velocity : velocityModel)
    {
        if (velocity <= 0)
        {
            throw std::invalid_argument("All velocities must be positive");
        }
    }
    for (int i = 0; i < static_cast<int> (velocityModel.size()) - 1; ++i)
    {
        if (velocityModel[i + 1] <= velocityModel[i])
        {
            throw std::invalid_argument("Velocity inversions not yet done");
        }
    }
    if (!std::is_sorted(interfaces.begin(), interfaces.end(),
                        [=](const double lhs, const double rhs)
                        {
                            return lhs < rhs;
                        }))
    {
        throw std::invalid_argument("Interfaces must increase with depth");
    }
    pImpl->mSlownessModel.resize(velocityModel.size());
    std::transform(velocityModel.begin(), velocityModel.end(),
                   pImpl->mSlownessModel.begin(),
                   [=](const double velocity)
                   {
                       return 1./velocity;
                   });
    pImpl->mInterfaces = interfaces;
}

bool LayerSolver::haveVelocityModel() const noexcept
{
    return !pImpl->mSlownessModel.empty();
}

/// Source depth
void LayerSolver::setSourceDepth(const double depth)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (depth < pImpl->mInterfaces[0])
    {
        throw std::invalid_argument("Source is in the air");
    }
    pImpl->mSourceDepth = depth; 
    pImpl->mSourceLayer = pImpl->getLayer(depth);
    pImpl->mHaveSourceDepth = true;
}

double LayerSolver::getSourceDepth() const
{
    if (!haveSourceDepth()){throw std::runtime_error("Source depth not set");}
    return pImpl->mSourceDepth;
}

bool LayerSolver::haveSourceDepth() const noexcept
{
    return pImpl->mHaveSourceDepth;
}

/// Station
void LayerSolver::setStationOffset(const double offset)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (offset < 0){throw std::invalid_argument("Offset must be non-negative");}
    pImpl->mStationOffset = offset;
    pImpl->mStationDepth = pImpl->mInterfaces.at(0);
    pImpl->mStationLayer = pImpl->getLayer(pImpl->mStationDepth);
    pImpl->mHaveStationOffsetAndDepth = true;
}

void LayerSolver::setStationOffsetAndDepth(const double offset,
                                           const double depth)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (depth < pImpl->mInterfaces[0])
    {   
        throw std::invalid_argument("Station is in the air");
    }
    if (offset < 0){throw std::invalid_argument("Offset must be non-negative");}
    pImpl->mStationOffset = offset;
    pImpl->mStationDepth = depth; 
    pImpl->mStationLayer = pImpl->getLayer(depth);
    pImpl->mHaveStationOffsetAndDepth = true;
#ifndef NDEBUG
    assert(pImpl->mStationLayer >= 0 &&
           pImpl->mStationLayer < static_cast<int> (pImpl->mInterfaces.size()));
#endif
}

bool LayerSolver::haveStationOffsetAndDepth() const noexcept
{
    return pImpl->mHaveStationOffsetAndDepth;
}

/// Get the ray paths
std::vector<Path2D> LayerSolver::getRayPaths() const
{
    if (!haveRayPaths()){throw std::runtime_error("Ray paths not computed");}
    return pImpl->mRayPaths;
}

const std::vector<Path2D> &LayerSolver::getRayPathsReference() const
{
    return *&pImpl->mRayPaths;
}

/// Have ray path?
bool LayerSolver::haveRayPaths() const noexcept
{
    return pImpl->mHaveRayPaths;
}

/// Solve
void LayerSolver::solve()
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (!haveSourceDepth()){throw std::runtime_error("Source depth not set");}
    if (!haveStationOffsetAndDepth())
    {
        throw std::runtime_error("Station offset not set");
    }
    Path2D rayPath;
    pImpl->mHaveRayPaths = false;
    pImpl->mRayPaths.clear();
    // Special case of just one layer
    if (pImpl->mSlownessModel.size() == 1)
    {
        pImpl->mRayPaths.push_back(pImpl->computeDirectRayInSameLayer());
        pImpl->mHaveRayPaths = true;
        return;
    }
    // Special case station directly above/below station
    if (pImpl->mStationOffset < std::numeric_limits<double>::epsilon())
    {
        int direction =+1;
        int sourceLayer = pImpl->mSourceLayer;
        int stationLayer = pImpl->mStationLayer;
        Point2D startPoint{0, pImpl->mSourceDepth};
        Point2D endPoint{pImpl->mStationOffset, pImpl->mStationDepth};
        // Same layer
        if (sourceLayer == stationLayer)
        {
            Segment2D segment;
            segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
            segment.setSlowness(pImpl->mSlownessModel.at(sourceLayer));
            segment.setVelocityModelCellIndex(sourceLayer);
            rayPath.open();
            rayPath.append(std::move(segment));
            rayPath.close();
            pImpl->mRayPaths.push_back(std::move(rayPath));
            pImpl->mHaveRayPaths = true;
            return;
        }
        // Trace `down' 
        int nextLayerIncrement = 1;
        if (pImpl->mSourceLayer > pImpl->mStationLayer)
        {
            direction =-1;
            nextLayerIncrement = 0;
        }
        auto currentLayer = sourceLayer;
        Point2D point0{startPoint};
        rayPath.open();
        while (currentLayer != stationLayer)
        {
            Segment2D segment;
            auto nextLayer = currentLayer + nextLayerIncrement;
            Point2D point1;
            point1.setPositionInX(0);
            point1.setPositionInZ(pImpl->mInterfaces.at(nextLayer));
            segment.setStartAndEndPoint(std::pair {point0, point1});
            segment.setVelocityModelCellIndex(currentLayer);
            segment.setSlowness(pImpl->mSlownessModel.at(currentLayer));
            if (segment.getLength() > 1.e-3)
            {
                rayPath.append(std::move(segment));
            }
            // Update `start point'
            point0 = point1;
            // Update layer
            currentLayer = currentLayer + direction;
            //std::cout << "currentLayer " << currentLayer << " " << sourceLayer << " " << stationLayer << std::endl;
        }
        // Last case - trace from interface
        Segment2D lastSegment;
        lastSegment.setStartAndEndPoint(std::pair {point0, endPoint});
        lastSegment.setVelocityModelCellIndex(pImpl->mStationLayer);
        lastSegment.setSlowness(pImpl->mSlownessModel.at(pImpl->mStationLayer));
        rayPath.append(std::move(lastSegment));
        rayPath.close();
        pImpl->mRayPaths.push_back(std::move(rayPath));
        pImpl->mHaveRayPaths = true;
        return;
    }
    // 
    // Solve general case
    std::vector<double> takeOffAngles;
    computeTakeOffAngles(45, 1.e-7, 90 - 1.e-7, &takeOffAngles);
    for (const auto &theta : takeOffAngles)
    {
        pImpl->shoot(theta);
//break;
    }
}
