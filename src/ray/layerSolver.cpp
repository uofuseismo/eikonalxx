#include <algorithm>
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

class LayerSolver::LayerSolverImpl
{
public:
    Path2D mRayPath;
    std::vector<double> mSlownessModel;
    std::vector<double> mInterfaces;
    double mStationDepth{0};
    double mSourceDepth{0};
    double mStationOffset{0};
    int mSourceLayer{0};
    int mStationLayer{0};
    bool mHaveStationOffsetAndDepth{false};
    bool mHaveSourceDepth{false};
    bool mHaveRayPath{false};
};

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
    if (depth >= pImpl->mInterfaces.back())
    {
        pImpl->mSourceLayer = static_cast<int> (pImpl->mInterfaces.size()) - 1;
    }
    else
    {
        pImpl->mSourceLayer
            = std::distance(
                  pImpl->mInterfaces.begin(),
                  std::lower_bound(pImpl->mInterfaces.begin(),
                                   pImpl->mInterfaces.end(),
                                   pImpl->mSourceDepth));
    }
    pImpl->mHaveSourceDepth = true;
#ifndef NDEBUG
    assert(pImpl->mSourceLayer >= 0 &&
           pImpl->mSourceLayer < static_cast<int> (pImpl->mInterfaces.size()));
#endif
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
    pImpl->mStationLayer = 0;
    pImpl->mHaveStationOffsetAndDepth = true;
#ifndef NDEBUG
    assert(pImpl->mStationLayer >= 0 &&
           pImpl->mStationLayer < static_cast<int> (pImpl->mInterfaces.size()));
#endif
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
    if (depth >= pImpl->mInterfaces.back())
    {
        pImpl->mStationLayer = static_cast<int> (pImpl->mInterfaces.size()) - 1;
    }
    else
    {
        pImpl->mStationLayer
            = std::distance(
                  pImpl->mInterfaces.begin(),
                  std::lower_bound(pImpl->mInterfaces.begin(),
                                   pImpl->mInterfaces.end(),
                                   pImpl->mSourceDepth));
    }
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

/// Have ray path?
bool LayerSolver::haveRayPath() const noexcept
{
    return pImpl->mHaveRayPath;
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
    pImpl->mHaveRayPath = false;
    pImpl->mRayPath.clear();
    // Special case
    if (pImpl->mSlownessModel.size() == 1)
    {
        Segment2D segment;
        Point2D startPoint{0, pImpl->mSourceDepth};
        Point2D endPoint{pImpl->mStationOffset, pImpl->mStationDepth};
        segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
        segment.setSlowness(pImpl->mSlownessModel[0]);
        segment.setVelocityModelCellIndex(0);
        pImpl->mRayPath.open();
        pImpl->mRayPath.append(std::move(segment));
        pImpl->mRayPath.close();
        pImpl->mHaveRayPath = true;
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
            segment.setSlowness(pImpl->mSlownessModel[sourceLayer]);
            segment.setVelocityModelCellIndex(sourceLayer);
            pImpl->mRayPath.open();
            pImpl->mRayPath.append(std::move(segment));
            pImpl->mRayPath.close();
            pImpl->mHaveRayPath = true;
            return;
        }
        // Trace `down' 
        if (pImpl->mSourceLayer > pImpl->mStationLayer)
        {
            direction =-1;
        }
        bool isFirst{true};
        auto currentLayer = sourceLayer;
        Segment2D segment;
        pImpl->mRayPath.open();
        Point2D point0{startPoint};
        Point2D point1;
        while (currentLayer != stationLayer)
        {
            if (direction ==+1)
            {
                if (!isFirst)
                {
                    point0.setPositionInX(0);
                    point0.setPositionInZ(pImpl->mInterfaces[currentLayer + 1]);
                }
                else
                {
                    isFirst = false;
                }
                point1.setPositionInX(0);
                point1.setPositionInZ(pImpl->mInterfaces[currentLayer]);
            }
            else
            {
                if (!isFirst)
                {   
                    point0.setPositionInX(0);
                    point0.setPositionInZ(pImpl->mInterfaces[currentLayer]);
                }
                else
                {   
                    isFirst = false;
                }   
                point1.setPositionInX(0);
                point1.setPositionInZ(pImpl->mInterfaces[currentLayer + 1]);
            } 
            segment.setStartAndEndPoint(std::pair {point0, point1});
            segment.setSlowness(pImpl->mSlownessModel[currentLayer]); 
            pImpl->mRayPath.append(segment);
            // Update `start point'
            point0 = point1;
        }
        // Last case - trace from interface
        point1 = endPoint;
        segment.setStartAndEndPoint(std::pair {point0, point1});
        segment.setSlowness(pImpl->mSlownessModel[pImpl->mStationLayer]);
        pImpl->mRayPath.append(segment);
        pImpl->mRayPath.close();
        pImpl->mHaveRayPath = true;
        return;
    }
}
