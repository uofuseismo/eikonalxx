#include <iostream>
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
    Overshot,     /*!< The ray exited the medium after the station offset. */
    ExitedBottom, /*!< The ray exited the bottom of the medium. */
    ExitedTop,    /*!< The ray exited the top of the medium. */
};

struct Path
{

};

}

class LayerSolver::LayerSolverImpl
{
public:
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
    void shoot(const double takeOffAngle)
    {
        constexpr int maxSegments{std::numeric_limits<int>::max()};
        constexpr double degreesToRadians{M_PI/180};
        auto currentAngle = takeOffAngle;
        auto currentLayer = mSourceLayer;
        for (int kSegment = 0; kSegment < maxSegments; ++kSegment)
        {
            int nextLayer =+1;
            if (currentAngle > 90)
            {
                nextLayer =-1;
            }
            // 
            // Ray exited top of medium.
            if (currentLayer + nextLayer < 0)
            {
                break;                
            }
        }
    }
    std::vector<Path2D> mRayPaths;
    std::vector<double> mSlownessModel;
    std::vector<double> mInterfaces;
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
    // Special case
    if (pImpl->mSlownessModel.size() == 1)
    {
        Segment2D segment;
        Point2D startPoint{0, pImpl->mSourceDepth};
        Point2D endPoint{pImpl->mStationOffset, pImpl->mStationDepth};
        segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
        segment.setSlowness(pImpl->mSlownessModel[0]);
        segment.setVelocityModelCellIndex(0);
        rayPath.open();
        rayPath.append(std::move(segment));
        rayPath.close();
        pImpl->mRayPaths.push_back(std::move(rayPath));
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
}
