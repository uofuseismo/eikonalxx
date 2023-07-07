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
#include "layerTracer.hpp"

/// 1.e-4 changes about 1 m for every 500 km
#define MIN_DOWNGOING_ANGLE 1.e-4
#define MAX_DOWNGOING_ANGLE 89.9999


using namespace EikonalXX::Ray;

namespace
{

enum class RayType
{
    Direct,
    Reflected,
    Refracted
};

//enum class ReturnCode
//{
//    Undershot,    /*!< The ray exited the of the medium before the station offset. */
//    Overshot,     /*!< The ray exited the medium after the station offset. */
//    ExitedBottom, /*!< The ray exited the bottom of the medium. */
//    ExitedRight,  /*!< The ray exited the right side of the medium. */
//    Hit,          /*!< The ray hit its intended target to within some tolerance. */
//    RayDoesNotTurn,
//    RayTurnsTooEarly
//};

/// Traces down then up where the source/receiver are at the same depth
/// ------------- Interface 0
/// *          x    Slowness 0       
/// ------------- Interface 1
///   \      /      Slowness 1
///-------------- Interface 2
///     \__/        Slowness 2
///-------------- Interace 3
///                 Slowness 3

/// Traces up:
/// ------------------------
/// -                   x  -
/// ------------------------
/// -                 /    -
/// -                /     -
/// ------------------------
/// -              /       -
/// -            /         -
/// -          *           -
/// ------------------------
ReturnCode traceDirectUp(const int stationLayer,
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
    path->clear();
    if (takeOffAngle >= 90)
    {
        throw std::invalid_argument("Ray must be going up");
    }
    // Turn ray around and make it go up
    auto takeOffAngleRadians = (180 - takeOffAngle)*(M_PI/180); 
    // Trace to top layer
    auto nLayers = static_cast<int> (interfaces.size());
    double x0{0};
    double z0{sourceDepth};
    double z1 = interfaces[startIndex]; 
    double dz = sourceDepth - z1;
#ifndef NDEBUG
    assert(dz >= 0); 
#endif
    double x1 = 0;
    if (dz > 0)
    {
        x1 = dz*std::tan(takeOffAngleRadians);
    }
    // Trace out first layer
    std::vector<::Segment> segments;
    segments.reserve(2*nLayers);
    ::Segment firstSegment{x0, z0, x1, z1, 
                           slownesses.at(startIndex), startIndex};
    segments.push_back(std::move(firstSegment));
    for (int layer = startIndex; layer >= endIndex + 1; --layer)
    {
         
    }
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
    /// For when the source and station are in the same layer this computes
    /// the direct ray. 
    [[nodiscard]] Path2D computeDirectRaySameLayer() const
    {
#ifndef NDEBUG
        assert(mSourceLayer == mStationLayer);
#endif
        Path2D path;
#ifndef NDEBUG
        auto returnCode =
#endif
        ::traceDirectSameLayer(mAugmentedInterfaces,
                               mAugmentedSlownesses,
                               mSourceLayer,
                               mSourceDepth,
                               mStationLayer,
                               mStationOffset,
                               mStationDepth,
                               &path);
#ifndef NDEBUG 
        assert(returnCode == ReturnCode::Hit);
#endif
        return path;
    }
    /// For the whole or half space problem
    [[nodiscard]] Path2D computeWholeSpace() const
    {
#ifndef NDEBUG
        assert(mSourceLayer == 0);
        assert(mStationLayer == mSourceLayer);
#endif
        Path2D rayPath;
#ifndef NDEBUG
        auto returnCode =
#endif
        ::traceWholeSpace(mSlownessModel[mSourceLayer],
                          mSourceDepth,
                          mStationOffset,
                          mStationDepth,
                          &rayPath);
#ifndef NDEBUG
        assert(returnCode == ReturnCode::Hit);
#endif
        return rayPath;
    }
    /// For when the source and station are vertically aligned this computes
    /// the direct ray (assuming the source is below the station) and the
    /// reflected ray paths.
    [[nodiscard]] std::vector<Path2D> computeVerticalRayPaths() const
    {
        std::vector<Path2D> rayPaths;
        // Trace the direct ray path upwards 
        if (mSourceDepth >= mStationDepth)
        {
            constexpr double upTakeOffAngle{180};
            std::vector<::Segment> segments;
#ifndef NDEBUG
            auto returnCode =
#endif
            ::traceDirect(mAugmentedInterfaces,
                          mAugmentedSlownesses,
                          upTakeOffAngle,
                          mSourceLayer,
                          mSourceDepth,
                          mStationLayer,
                          mStationOffset,
                          mStationDepth,
                          &segments,
                          mRayHitTolerance);
#ifndef NDEBUG 
            assert(returnCode == ReturnCode::Hit ||
                  returnCode == ReturnCode::UnderShot);
#endif
            if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
        }
        else // Same layer but station below
        {
            if (mSourceLayer == mStationLayer)
            {
                rayPaths.push_back(computeDirectRaySameLayer());
            }
        }
        // Loop through stack and bounce rays
        auto nLayers = static_cast<int> (mInterfaces.size());
        for (int layer = mSourceLayer; layer < nLayers - 1; ++layer)
        {
            std::vector<::Segment> segments;
#ifndef NDEBUG
            auto returnCode = 
#endif
                 ::traceVerticalReflectionDown(mAugmentedInterfaces,
                                               mAugmentedSlownesses,
                                               mSourceLayer,
                                               layer,
                                               mSourceDepth,
                                               mStationLayer,
                                               mStationDepth,
                                               mStationOffset,
                                               &segments,
                                               mRayHitTolerance);
#ifndef NDEBUG
            assert(returnCode == ReturnCode::Hit ||
                  returnCode == ReturnCode::UnderShot);
#endif
            if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
        }
        std::sort(rayPaths.begin(), rayPaths.end(),
                  [](const Path2D &lhs, const Path2D &rhs)
                  {
                      return lhs.getTravelTime() < rhs.getTravelTime();
                  });
        return rayPaths;
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
    /// @brief Shoots a ray with the given take-off angle.
    /// @note This cannot handle velocity inversions. 
    std::vector<Path2D> shoot(const double takeOffAngle,
                              const bool keepOnlyHits = true)
    {
        std::vector<Path2D> result;
        // Edge case - straight down or up
        if (takeOffAngle < 1.e-4 || std::abs(180 - takeOffAngle) < 1.e-4)
        {
            if (keepOnlyHits && mStationOffset > 1.e-4)
            {
                auto temporaryRays = computeVerticalRayPaths();
                result.reserve(temporaryRays.size());
                for (auto &ray : temporaryRays)
                {
                    auto nSegments = ray.size();
                    if (nSegments > 0)
                    {
                        auto offset
                             = ray.at(nSegments - 1).getEndPoint().getPositionInX();
                        if (std::abs(offset - mStationOffset) < mRayHitTolerance)
                        {
                            result.push_back(std::move(ray));
                        }
                    }
                }
            }
            else
            {
                result = computeVerticalRayPaths();
            } 
            return result;
        }
        // 90 degree take-off angles won't converge
        if (std::abs(takeOffAngle - 90) < 1.e-4){return result;}
        // Up-going rays
        if (takeOffAngle > 90)
        {
            // Not bouncing off above layers so just quit now
            if (mStationDepth < mSourceDepth){return result;}
            std::vector<::Segment> segments;
            try
            {
                auto returnCode = ::traceDirect(mAugmentedInterfaces,
                                                mAugmentedSlownesses,
                                                takeOffAngle,
                                                mSourceLayer,
                                                mSourceDepth,
                                                mStationLayer,
                                                mStationOffset,
                                                mStationDepth,
                                                &segments,
                                                mRayHitTolerance); 
                if (returnCode == ReturnCode::Hit ||
                    (!keepOnlyHits && returnCode == ReturnCode::UnderShot) ||
                    (!keepOnlyHits && returnCode == ReturnCode::OverShot))
                {
                    auto rayPath = ::toRayPath(segments);
                    //std::cout << "Takeoff angle, travel time: " << takeOffAngle << "," << rayPath.getTravelTime() << std::endl;
                    result.push_back(rayPath);
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << std::endl;
            }
        }
        else
        {
            // Trace down through the velocity model stack and tabulate the
            // ray paths
            auto nLayers = static_cast<int> (mInterfaces.size());
            for (int endLayer = mSourceLayer; endLayer < nLayers - 1; ++endLayer)
            {
                std::vector<::Segment> segments;
                try
                {
                    auto returnCode = ::traceDownThenUp(mAugmentedInterfaces,
                                                        mAugmentedSlownesses,
                                                        takeOffAngle,
                                                        mSourceLayer,
                                                        mSourceDepth,
                                                        mStationLayer,
                                                        mStationOffset,
                                                        mStationDepth,
                                                        endLayer,
                                                        &segments,
                                                        mRayHitTolerance);
                    if (returnCode == ReturnCode::Hit ||
                       (!keepOnlyHits && returnCode == ReturnCode::UnderShot) ||
                       (!keepOnlyHits && returnCode == ReturnCode::OverShot))
                    {
                        auto rayPath = ::toRayPath(segments);
/*
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
                        std::cout << "Takeoff angle, travel time: " << takeOffAngle << "," << rayPath.getTravelTime() << std::endl;
*/
                        result.push_back(rayPath);
                    }
                }
                catch (const std::exception &e)
                {
                    std::cerr << e.what() << std::endl;
                }
            } // Loop on stack of layers
        } // End check on upgoing vs downgoing
        // Sort the rays in increasing order based on the travel times
        if (result.size() > 1)
        {
            std::sort(result.begin(), result.end(), 
                      [](const Path2D &lhs, const Path2D &rhs)
                      {
                         return lhs.getTravelTime() < rhs.getTravelTime();
                      });
        }
        return result;
    }
    std::pair<double, Path2D>
         quadraticFitSearch(const std::pair<double, double> &fa,
                            const std::pair<double, double> &fb,
                            const std::pair<double, double> &fc,
                            const int maxFunctionEvaluations = 5)
    {
        constexpr bool keepOnlyHits{true};
        auto a  = fa.first;
        auto ya = fa.second;        
        auto b  = fb.first;
        auto yb = fb.second;
        auto c  = fc.first;
        auto yc = fc.second;
        auto xStar = 1./3.*(a + b + c);
        auto fStar = std::max(ya, std::max(yb, yc));
        Path2D fPath;
        for (int k = 0; k < maxFunctionEvaluations; ++k)
        {
#ifndef NDEBUG
            assert(a < b);
            assert(b < c);
#endif
            auto a2 = a*a;
            auto b2 = b*b;
            auto c2 = c*c;
            auto xNum = ya*(b2 - c2) + yb*(c2 - a2) + yc*(a2 - b2);
            auto xDen = 2*(ya*(b - c) + yb*(c - a) +  yc*(a - b));
            xStar = 1./3.*(a + b + c);
            if (std::abs(xDen) > 1.e-14)
            {
                xStar = xNum/xDen;
            }
            else
            {
                std::cerr << "Division by 0: " << xNum << ","
                          << xDen << std::endl;
            }
            auto paths = shoot(xStar, keepOnlyHits);
            if (paths.empty())
            {
                throw std::runtime_error("No paths converged");
            }
            fPath = paths.at(0);
            fStar = fPath.getTravelTime(); 
            if (xStar > b)
            {
                if (fStar > yb)
                {
                    c = xStar;
                    yc = fStar; 
                }
                else
                {
                    a  = b;
                    ya = yb;
                    b  = xStar;
                    yb = fStar;
                }
            }
            else if (xStar < b)
            {
                if (fStar > yb)
                {
                    a = xStar;
                    ya = fStar;
                }
                else
                {
                    c  = b;
                    yc = yb;
                    b  = xStar;
                    yb = fStar;
                }
            }
            else
            {
                return std::pair {b, fPath};
            }
        }
        return std::pair {xStar, fPath};
    } 
    [[nodiscard]] Path2D optimizeDownGoing()
    {
        Path2D optimalPath;
        // Initial bracketing
        int nInitialAngles{89};
        std::vector<double> takeOffAngles;
        computeTakeOffAngles(nInitialAngles,
                             MIN_DOWNGOING_ANGLE, MAX_DOWNGOING_ANGLE,
                             &takeOffAngles);
        auto dAngle2 = std::abs(takeOffAngles.at(1) - takeOffAngles.at(0))/2;
        std::vector<std::pair<double, double>> angleTravelTime;
        for (auto &takeOffAngle : takeOffAngles)
        {
            auto work = shoot(takeOffAngle);
            if (!work.empty())
            {
                angleTravelTime.push_back(std::pair {takeOffAngle,
                                                     work[0].getTravelTime()}); 
            }
        }
        // Probably too close
        if (angleTravelTime.empty())
        {
            //std::cout << "No downgoing ray paths converged" << std::endl;
            return optimalPath;
        }
        // Identify each candidate solution region
        std::vector<std::pair<int, int>> brackets;
        if (angleTravelTime.size() == 1)
        {
            auto theta0 = std::max(MIN_DOWNGOING_ANGLE,
                                   angleTravelTime[0].first - dAngle2);
            auto theta1 = std::min(MAX_DOWNGOING_ANGLE,
                                   angleTravelTime[0].first + dAngle2);
            brackets.push_back(std::pair {theta0, theta1});
        }
        else
        {
            auto nAngles = static_cast<int> (angleTravelTime.size());
            for (int i = 1; i < nAngles - 1; ++i)
            {
                auto thisTravelTime = angleTravelTime[i].second;
/*
                if (i == 0)
                {
                    if (angleTravelTime[i + 1].second > thisTravelTime)
                    {
                        auto theta0
                            = std::max(MIN_DOWNGOING_ANGLE,
                                       angleTravelTime[i].first - dAngle2);
                        auto theta1 = angleTravelTime[i + 1].first;
                        brackets.push_back(std::pair {theta0, theta1});
                    }
                }
*/
                if (i > 0 && i < nAngles - 1) 
                {
                    if (angleTravelTime[i - 1].second > thisTravelTime &&
                        angleTravelTime[i + 1].second > thisTravelTime)
                    {
                        brackets.push_back(std::pair {i - 1, i + 1});
                    }
                }
/*
                else if (i == nAngles - 1)
                {
                    if (angleTravelTime[i - 1].second > thisTravelTime)
                    {
                        auto theta0 = angleTravelTime[i - 1].first;
                        auto theta1
                            = std::min(MAX_DOWNGOING_ANGLE,
                                       angleTravelTime[i].first + dAngle2);
                        brackets.push_back(std::pair {theta0, theta1});
                    }
                }
*/
#ifndef NDEBUG
                else
                {
                    assert(false);
                }
#endif
            }
        }
        // Search through each bracket
        double optimumTakeOffAngle{0};
        double minimumTravelTime{std::numeric_limits<double>::max()}; 
        bool failed{false};
        for (const auto &bracket : brackets)
        {
             auto ia = bracket.first;
             auto ib = ia + 1;
             auto ic = bracket.second;
             try
             {
                 auto [takeOffAngle, path]
                     = quadraticFitSearch(angleTravelTime.at(ia),
                                          angleTravelTime.at(ib),
                                          angleTravelTime.at(ic),
                                          5);
                 //std::cout << "Take-off angle, travel time: " << takeOffAngle << "," << path.getTravelTime() << std::endl;
                 if (path.getTravelTime() < minimumTravelTime)
                 {
                     optimalPath = path;
                     optimumTakeOffAngle = takeOffAngle;
                     minimumTravelTime = path.getTravelTime();
                 }
             }
             catch (const std::exception &e)
             {
                 std::cerr << e.what() << std::endl;
             }
        }
        return optimalPath;
    }
    std::vector<Path2D> mRayPaths;
    std::vector<double> mSlownessModel;
    std::vector<double> mInterfaces;
    std::vector<double> mAugmentedSlownesses;
    std::vector<double> mAugmentedInterfaces;
    double mRayHitTolerance{1};
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
    pImpl->mSlownessModel = ::toSlownessVector(velocityModel);
    pImpl->mInterfaces = interfaces;
    pImpl->mAugmentedSlownesses
         = ::toSlownessVector(::augmentVelocityVector(velocityModel));
    pImpl->mAugmentedInterfaces
         = ::augmentInterfacesVector(pImpl->mInterfaces);
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

/// Shoot
std::vector<Path2D> LayerSolver::shoot(const double takeOffAngle,
                                       const bool keepOnlyHits)
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
    if (takeOffAngle < 0 || takeOffAngle > 180)
    {
        throw std::invalid_argument("Take-off angle must be in range [0,180]");
    }
    return pImpl->shoot(takeOffAngle, keepOnlyHits);
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
        pImpl->mRayPaths.push_back(pImpl->computeWholeSpace());
        pImpl->mHaveRayPaths = true;
        return;
    }
    // Special case station directly above/below station
    if (pImpl->mStationOffset < std::numeric_limits<double>::epsilon())
    {
        pImpl->mRayPaths = pImpl->computeVerticalRayPaths();
        pImpl->mHaveRayPaths = true;
        return;
    }
    // General case where we have to optimize
//pImpl->shoot(39);//5);
//getchar();
    // Are we in the same layer?  Then just shoot direct
    if (pImpl->mSourceLayer == pImpl->mStationLayer)
    {
        pImpl->mRayPaths.push_back(pImpl->computeDirectRaySameLayer());
    }
    // Solve general case
    auto optimalDownGoingPath = pImpl->optimizeDownGoing();
    if (optimalDownGoingPath.size() > 0)
    {
/*
        for (const auto &p : optimalDownGoingPath)
        {
            std::cout << std::setprecision(12) << p.getStartPoint().getPositionInX() << ","
                      << p.getStartPoint().getPositionInZ() << ","
                      << p.getEndPoint().getPositionInX() << ","
                      << p.getEndPoint().getPositionInZ() << ","
                      << p.getVelocity() << std::endl;
        }
*/
        pImpl->mRayPaths.push_back(std::move(optimalDownGoingPath));
    }
    // Sort these
    std::sort(pImpl->mRayPaths.begin(), pImpl->mRayPaths.end(), 
              [&](const Path2D &lhs, const Path2D &rhs)
              {
                  return lhs.getTravelTime() < rhs.getTravelTime();
              });
    pImpl->mHaveRayPaths = true;
/*
    std::vector<double> takeOffAngles;
    computeTakeOffAngles(45, 1.e-7, 90 - 1.e-7, &takeOffAngles);
    for (const auto &theta : takeOffAngles)
    {
        pImpl->shoot(theta);
//break;
    }
*/
}
