#include <iomanip>
#include <iostream>
#include <limits>
#include <CL/sycl.hpp>
#include <algorithm>
#include <map>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#if __has_include(<pstl/execution>)
   #include <pstl/execution>
   #include <pstl/algorithm>
   #define USE_PSTL
#endif
//#if __has_include(<oneapi/dpl/execution>)
//   #include <oneapi/dpl/execution>
//   #include <oneapi/dpl/algorithm>
//   #define ONE_API
//#else
//#endif
#include "eikonalxx/ray/gradientTracer2d.hpp"
#include "eikonalxx/ray/gradientTracerOptions.hpp"
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/station2d.hpp"
#include "eikonalxx/io/vtkLines2d.hpp"
#include "private/grid.hpp"
#include "bilinear.hpp"

using namespace EikonalXX::Ray;

namespace
{

struct Segment
{
    double x0;
    double z0;
    double x1;
    double z1;
    int iCellX0;
    int iCellZ0;
    int iCellX1;
    int iCellZ1;
};

void reverseSegments(std::vector<::Segment> &segments)
{
    std::reverse(segments.begin(), segments.end());
    for (auto &segment : segments)
    {
        std::swap(segment.x0,      segment.x1);
        std::swap(segment.z0,      segment.z1);
        std::swap(segment.iCellX0, segment.iCellX1);
        std::swap(segment.iCellZ0, segment.iCellZ1);
    }
}

[[nodiscard]] bool geometryMatches(
    const EikonalXX::Geometry2D &lhs,
    const EikonalXX::Geometry2D &rhs)
{
    if (lhs.getNumberOfGridPointsInX() != rhs.getNumberOfGridPointsInX())
    {
        return false;
    }
    if (lhs.getNumberOfGridPointsInZ() != rhs.getNumberOfGridPointsInZ())
    {
        return false;
    }
    // Grid spacing should be good to millimeters
    if (std::abs(lhs.getGridSpacingInX() - rhs.getGridSpacingInX()) > 1.e-3)
    {
        return false;
    }
    if (std::abs(lhs.getGridSpacingInZ() - rhs.getGridSpacingInZ()) > 1.e-3)
    {
        return false;
    }
    // Origin should be good to centimeters
    if (std::abs(lhs.getOriginInX() - rhs.getOriginInX()) > 1.e-2)
    {
        return false;
    }
    if (std::abs(lhs.getOriginInZ() - rhs.getOriginInZ()) > 1.e-2)
    {
        return false;
    }

    return true;
}

/*
template<typename T>
T bilinear(const T x, const T z,
           const T x1, const T x2,
           const T z1, const T z2,
           const T f00, const T f01,
           const T f10, const T f11,
           const T dxi, const T dzi)
{
    T x2mx = x2 - x;
    T xmx1 = x - x1;
    T fxz1 = dxi*(x2mx*f00 + xmx1*f10);
    T fxz2 = dxi*(x2mx*f01 + xmx1*f11);

    T z2mz = z2 - z;
    T zmz1 = z - z1; 
    return dzi*(z2mz*fxz1 + zmz1*fxz2);
}
*/

/// @result The point where this line intersects the next cell wall.
/// The idea for intersection is pretty straightforward.  We have a 
/// point and a gradient, hence we are moving along the linear 
/// trajectory
///   \textbf{f}(t) = \textbf{o} + \textbf{g} t 
/// where \textbf{o} is the origin of the vector, \textbf{g} is 
/// direction given by the gradient field, and t is some unknown value.
/// Now, we find where \textbf{v} intersects the edge of the cell.  
/// Really, all we are now doing is finding the point of intersection
/// with some existence conditions:
///   Left:   o_x + g_x t_l = x_0 -> t_l = (x_0 - o_x)/g_x; g_x < 0
///   Right:  o_x + g_x t_r = x_1 -> t_r = (x_1 - o_x)/g_x; g_x > 0
///   Bottom: o_z + g_z t_b = z_0 -> t_b = (z_0 - o_z)/g_z; g_z < 0
///   Top:    o_z + g_z t_t = z_1 -> t_t = (z_1 - o_z)/g_z; g_z > 0
/// Barring strict equality of a gradient component with zero, we
/// will have two candidate solutions for t; one which intersects the
/// appropriate x plane and one which intersects the appropriate z plane.
/// The t we should pick should be the t that yields the shortest line
/// - i.e., intersects the closest plane.  Therefore, if 
///     | \textbf{o} + \textbf{g} t_i | < | \textbf{o} + \textbf{g} t_j |
/// we would select t_i, otherwise, we select t_j. Following this 
void intersect(const int iCellX0, const int iCellZ0,
               const double ox, const double oz,
               const double gx, const double gz,
               const double dx, const double dz,
               double *ex, double *ez)
{
    *ex = ox;
    *ez = oz;
    double ti{1.e10};
    double tj{1.e10};
    // Bracket the position with horizontal and vertical planes
    double x0 = iCellX0*dx;
    double x1 = (iCellX0 + 1)*dx;
    double z0 = iCellZ0*dz;
    double z1 = (iCellZ0 + 1)*dz;
#ifndef NDEBUG
//    assert(ox >= x0 && ox <= x1);
//    assert(oz >= z0 && oz <= z1);
#endif
    /// Existence
    if (gx > 0)
    {
        ti = (x1 - ox)/gx;
    }
    else if (gx < 0)
    {
        ti = (x0 - ox)/gx;
    }
    if (gz > 0)
    {
        tj = (z1 - oz)/gz;
    }
    else if (gz < 0)
    {
        tj = (z0 - oz)/gz;
    }
    // Make my lines
    double fxi = ox + gx*ti;
    double fzi = oz + gz*ti;
    double fxj = ox + gx*tj;
    double fzj = oz + gz*tj;
    // Distance squared is good enough
    double li = (fxi - ox)*(fxi - ox) + (fzi - oz)*(fzi - oz);
    double lj = (fxj - ox)*(fxj - ox) + (fzj - oz)*(fzj - oz);
    if (li < lj)
    {
        *ex = fxi;
        *ez = fzi;
    }
    else
    {
        *ex = fxj;
        *ez = fzj;
    }
}

}

class GradientTracer2D::GradientTracer2DImpl
{
public:
    EikonalXX::Geometry2D mGeometry;
    EikonalXX::Ray::GradientTracerOptions mOptions;
    std::vector<EikonalXX::Station2D> mStations;
    std::vector<EikonalXX::Ray::Path2D> mPaths;
    bool mInitialized{false};
    bool mHaveRayPaths{false};
};

/// Constructor
GradientTracer2D::GradientTracer2D() :
    pImpl(std::make_unique<GradientTracer2DImpl> ())
{
}

/// Copy constructor
GradientTracer2D::GradientTracer2D(const GradientTracer2D &tracer)
{
    *this = tracer;
}

/// Move constructor
GradientTracer2D::GradientTracer2D(GradientTracer2D &&tracer) noexcept
{
    *this = std::move(tracer);
}

/// Destructor
GradientTracer2D::~GradientTracer2D() = default;

/// Copy assignment
GradientTracer2D& GradientTracer2D::operator=(const GradientTracer2D &tracer)
{
    if (&tracer == this){return *this;}
    pImpl = std::make_unique<GradientTracer2DImpl> (*tracer.pImpl);
    return *this;
}

/// Move assignment
GradientTracer2D&
GradientTracer2D::operator=(GradientTracer2D &&tracer) noexcept
{
    if (&tracer == this){return *this;}
    pImpl = std::move(tracer.pImpl);
    return *this;
}

/// Set stations
void GradientTracer2D::setStations(
    const std::vector<EikonalXX::Station2D> &stations)
{
    if (!isInitialized()){throw std::invalid_argument("Class not initialized");}
    for (const auto &station : stations)
    {
        if (!station.haveGeometry())
        {
            throw std::invalid_argument("Station geometry not set");
        }
        if (!station.haveLocationInX())
        {
            throw std::invalid_argument("Station x position not set");
        }
        if (!station.haveLocationInZ())
        {
            throw std::invalid_argument("Station z position not set");
        }
        if (!::geometryMatches(pImpl->mGeometry, station.getGeometry()))
        {
            throw std::invalid_argument("Inconsistent geometry");
        }
    }
    pImpl->mStations = stations;
    pImpl->mPaths.clear();
    pImpl->mHaveRayPaths = false;
}

/// Have stations?
bool GradientTracer2D::haveStations() const noexcept
{
    return !pImpl->mStations.empty();
}

/// Set the geometry
void GradientTracer2D::initialize(
    const EikonalXX::Ray::GradientTracerOptions &options,
    const EikonalXX::Geometry2D &geometry)
{
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("x grid spacing not set");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("z grid spacing not set");
    }
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Grid points in x not set");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Grid points in z not set");
    }
    pImpl->mOptions = options;
    pImpl->mGeometry = geometry;
    pImpl->mInitialized = true;
    pImpl->mPaths.clear();
    pImpl->mHaveRayPaths = false;
}

/// Initialized?
bool GradientTracer2D::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Trace
template<typename T>
void GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<T> &solver)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!solver.haveVelocityModel())
    {
        throw std::invalid_argument("Velocity model not set on solver");
    }
    if (!solver.haveSource())
    {
        throw std::invalid_argument("Source not on solver");
    }
    if (!solver.haveTravelTimeGradientField())
    {
        throw std::invalid_argument("Travel time gradient field not set");
    }
    if (!haveStations())
    {
        throw std::runtime_error("No stations were set");
    }
    // Initialize
    pImpl->mHaveRayPaths = false;
    pImpl->mPaths.clear();
    pImpl->mPaths.resize(pImpl->mStations.size());
    // Grid properties
    int nGridX = pImpl->mGeometry.getNumberOfGridPointsInX();
    int nGridZ = pImpl->mGeometry.getNumberOfGridPointsInZ();
    //int nCellX = pImpl->mGeometry.getNumberOfCellsInX();
    //int nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
    double dx = pImpl->mGeometry.getGridSpacingInX();
    double dz = pImpl->mGeometry.getGridSpacingInZ();
    double dxi = 1./dx;
    double dzi = 1./dz;
    double x0 = pImpl->mGeometry.getOriginInX();
    double z0 = pImpl->mGeometry.getOriginInZ();
    // Step properties
    auto radiusScaleFactors = pImpl->mOptions.getRadiusScaleFactor();
    double defaultStepLength
         = radiusScaleFactors.back().second*std::min(dx, dz);
    // Source properties
    auto source = solver.getSource();
    auto iSourceCellX = source.getCellInX();
    auto iSourceCellZ = source.getCellInZ();
    auto iSourceCell = source.getCell();
    auto xs = source.getOffsetInX();
    auto zs = source.getOffsetInZ();
    auto sourceSlowness
        = static_cast<double> (solver.getSlowness(iSourceCellX, iSourceCellZ));
    // Set the source point (we'll always use this when tracing)
    Point2D sourcePoint{x0 + xs, z0 + zs}; 
    const T* gradientPtr = solver.getTravelTimeGradientFieldPointer();
    // Loop on receivers
    for (int iStation = 0; iStation < static_cast<int> (pImpl->mStations.size()); ++iStation)
    {
        auto iStationCellX = pImpl->mStations[iStation].getCellInX();
        auto iStationCellZ = pImpl->mStations[iStation].getCellInZ();
        // Define the station piont (we'll also always use this when tracing)
        auto xr = pImpl->mStations[iStation].getOffsetInX();
        auto zr = pImpl->mStations[iStation].getOffsetInZ();
        Point2D stationPoint {x0 + xr, z0 + zr};
        auto iStationCell = pImpl->mStations[iStation].getCell();
        auto stationSlowness
            = static_cast<double> (solver.getSlowness(iStationCellX,
                                                      iStationCellZ));
        // Deal with an algorithm breakdown where the source/station are in
        // the same cell
        if (iSourceCell == iStationCell)
        {
            Segment2D segment;
            segment.setStartAndEndPoint(std::pair {sourcePoint, stationPoint});
            segment.setSlowness(
                solver.getSlowness(iStationCellX, iStationCellZ));
            segment.setVelocityModelCellIndex(iStationCell); 
            Path2D rayPath;
            rayPath.open();
            rayPath.append(std::move(segment));
            rayPath.close();
            pImpl->mPaths[iStation] = std::move(rayPath);
            continue; 
        }
        // Business as usual - march down the gradient
        double xRayStart = xr;
        double zRayStart = zr;
        double xRayEnd = 0;
        double zRayEnd = 0;
        //int iCurrentCellX = iStationCellX;
        //int iCurrentCellZ = iStationCellZ;
        //int iNextCellX = 0;
        //int iNextCellZ = 0;
        int iCellX0 =-1;
        int iCellZ0 =-1;
        double gx000 = 1.e10;
        double gz000 = 1.e10;
        double gx100 = 1.e10;
        double gz100 = 1.e10;
        double gx010 = 1.e10;
        double gz010 = 1.e10;
        double gx110 = 1.e10;
        double gz110 = 1.e10;
        double stepLength0 = defaultStepLength;
        bool mConverged{false};
        const int maxSegments{2500};
        int nCellsVisited = 1;
        std::vector<::Segment> raySegments;
        raySegments.reserve(3*(nGridX + nGridZ));
        for (int iSegment = 0; iSegment < maxSegments; ++iSegment)
        {
            auto iCellXStart = static_cast<int> (xRayStart/dx);
            auto iCellZStart = static_cast<int> (zRayStart/dz);
            auto stepLength = stepLength0;
            auto gx00 = gx000;
            auto gz00 = gz000;
            auto gx10 = gx100;
            auto gz10 = gz100;
            auto gx01 = gx010;
            auto gz01 = gz010;
            auto gx11 = gx110;
            auto gz11 = gz110;
            // The ray is starting in a different cell than last time
            if (iCellXStart != iCellX0 || iCellZStart != iCellZ0)
            {
                // Find appropriate step length
                auto iCellDistance
                    = std::min(std::abs(iSourceCellX - iCellXStart),
                               std::abs(iSourceCellZ - iCellZStart));
                stepLength = radiusScaleFactors.back().second*std::min(dx, dz);
                for (const auto &rsf : radiusScaleFactors)
                {
                    if (rsf.first >= iCellDistance)
                    {
                        stepLength = rsf.second*std::min(dx, dz);
                        break;
                    }   
                }
                // Extract gradient and interpolate
                auto i00 = 2*::gridToIndex(nGridX, iCellXStart,     iCellZStart);
                auto i10 = 2*::gridToIndex(nGridX, iCellXStart + 1, iCellZStart);
                auto i11 = 2*::gridToIndex(nGridX, iCellXStart + 1, iCellZStart + 1);
                auto i01 = 2*::gridToIndex(nGridX, iCellXStart,     iCellZStart + 1);
                gx00 = static_cast<double> (gradientPtr[i00]);
                gz00 = static_cast<double> (gradientPtr[i00 + 1]);
                gx10 = static_cast<double> (gradientPtr[i10]);
                gz10 = static_cast<double> (gradientPtr[i10 + 1]);
                gx01 = static_cast<double> (gradientPtr[i01]);
                gz01 = static_cast<double> (gradientPtr[i01 + 1]);
                gx11 = static_cast<double> (gradientPtr[i11]);
                gz11 = static_cast<double> (gradientPtr[i11 + 1]);
                nCellsVisited = nCellsVisited + 1;
            }
            auto gx = ::bilinear(xRayStart, zRayStart,
                                 iCellXStart*dx, (iCellXStart + 1)*dx,
                                 iCellZStart*dz, (iCellZStart + 1)*dz,
                                 gx00, gx01,
                                 gx10, gx11,
                                 dxi, dzi);
            auto gz = ::bilinear(xRayStart, zRayStart,
                                 iCellXStart*dx, (iCellXStart + 1)*dx,
                                 iCellZStart*dz, (iCellZStart + 1)*dz,
                                 gz00, gz01,
                                 gz10, gz11,
                                 dxi, dzi);
            auto gNorm = std::hypot(gx, gz);
            auto gStep = stepLength/gNorm;
            // March a small step opposite direction of gradient
            xRayEnd = xRayStart - gx*gStep;
            zRayEnd = zRayStart - gz*gStep;
            auto iCellXEnd = static_cast<int> (xRayEnd/dx);
            auto iCellZEnd = static_cast<int> (zRayEnd/dz);
            // Update temporary segments
            ::Segment thisSegment{xRayStart, zRayStart,
                                  xRayEnd, zRayEnd,
                                  iCellXStart, iCellZStart,
                                  iCellXEnd,   iCellZEnd};
            raySegments.push_back(std::move(thisSegment));
            // Look ahead - did the ray finish its journal the source?
            if (std::abs(iCellXEnd - iSourceCellX) < 2 &&
                std::abs(iCellZEnd - iSourceCellZ) < 2)
            {
                ::Segment lastSegment{xRayEnd, zRayEnd,
                                      xs, zs,
                                      iCellXEnd, iCellZEnd,
                                      iSourceCellX, iSourceCellZ};
                raySegments.push_back(std::move(lastSegment));
                mConverged = true;
                break;
            }
            // Update
            stepLength0 = stepLength;
            xRayStart = xRayEnd;
            zRayStart = zRayEnd;
            iCellX0 = iCellXStart;
            iCellZ0 = iCellZStart;
            gx000 = gx00;
            gz000 = gz00;
            gx100 = gx10;
            gz100 = gz10;
            gx010 = gx01;
            gz010 = gz01;
            gx110 = gx11;
            gz110 = gz11;
        }
        //if (mConverged){std::cout << "Converged!" << std::endl;}
        if (!mConverged)
        {
            std::cerr << "Ray for station: "
                      << pImpl->mStations[iStation].getName()
                      << " did not converge to source" << std::endl;
            continue;
        }
        // Traced from receiver to source - flip that around
        ::reverseSegments(raySegments);
        // Start the process of healing the ray
        Path2D rayPath;
        rayPath.open();
        auto nSegments = static_cast<int> (raySegments.size());
        double xRay0 = raySegments[0].x0;
        double zRay0 = raySegments[0].z0;
        bool gotSource = false;
        for (int iSegment = 0; iSegment < nSegments; ++iSegment)
        {
            auto iCellX0 = raySegments[iSegment].iCellX0;
            auto iCellZ0 = raySegments[iSegment].iCellZ0;
            auto iCellX1 = raySegments[iSegment].iCellX1;
            auto iCellZ1 = raySegments[iSegment].iCellZ1;
            // We've hit the end
            if (iSegment == nSegments - 1)
            {
                Point2D startPoint{x0 + xRay0, z0 + zRay0};
                Segment2D segment;
                segment.setStartAndEndPoint(
                    std::pair{startPoint, stationPoint} );
                segment.setSlowness(stationSlowness);
                rayPath.append(std::move(segment));
                gotSource = true;
                break;
            }
            // We've crossed a cell boundary
            if (iCellX0 != iCellX1 || iCellZ0 != iCellZ1)
            {
                double gx = raySegments[iSegment].x1 - raySegments[iSegment].x0;
                double gz = raySegments[iSegment].z1 - raySegments[iSegment].z0;
                auto gnorm = std::hypot(gx, gz);
                gx = gx/gnorm;
                gz = gz/gnorm;
                double xRay1{0};
                double zRay1{0};
                //int iCellRayX0 = std::min(nCellX - 1, static_cast<int> (xRay0/dx));
                //int iCellRayZ0 = std::min(nCellZ - 1, static_cast<int> (zRay0/dz));
//std::cout << xRay0 << " " << zRay0 << std::endl;
                ::intersect(iCellX0, iCellZ0,
                            xRay0, zRay0,
                            gx, gz,
                            dx, dz,
                            &xRay1, &zRay1);
                Point2D startPoint{x0 + xRay0, z0 + zRay0};
                Point2D endPoint{x0 + xRay1, z0 + zRay1};
//std::cout << iSegment << " " << iCellX0 << " " << iCellZ0 << " " << x0 + xRay0 << " " << z0 + zRay0 << " " << x0 + xRay1 << " " << z0 + zRay1 << std::endl;
                //std::cout << xRay0 << " " << zRay0 << " " << xRay1 << " " << zRay1 << std::endl;
                Segment2D segment;
                segment.setStartAndEndPoint( std::pair{startPoint, endPoint} );
                if (iSegment == 0)
                {
                    segment.setSlowness(sourceSlowness);
                }
                else
                {
                    segment.setSlowness(solver.getSlowness(iCellX0, iCellZ0));
                }
                rayPath.append(std::move(segment));
                // Update
                xRay0 = xRay1;
                zRay0 = zRay1;
            }
        } // Loop
        if (!gotSource)
        {
            std::cerr << "Source point not in path" << std::endl;
        }
        rayPath.close();
        pImpl->mPaths[iStation] = std::move(rayPath);
//getchar();
        //getchar();
    } // Loop on stations
    pImpl->mHaveRayPaths = true;
}

bool GradientTracer2D::haveRayPaths() const noexcept
{
    return pImpl->mHaveRayPaths;
}

void GradientTracer2D::writeVTK(const std::string &fileName,
                                const std::string &title)
{
    if (!haveRayPaths()){throw std::runtime_error("Ray paths not computed");}
    EikonalXX::IO::VTKLines2D vtkWriter;
    vtkWriter.open(fileName, pImpl->mGeometry, title);
    if (!pImpl->mPaths.empty())
    {
        std::vector<std::vector<std::pair<double, double>>> lines;
        lines.reserve(pImpl->mPaths.size());
        for (int iPath = 0;
             iPath < static_cast<int> (pImpl->mPaths.size()); ++iPath)
        {
            auto path = pImpl->mPaths.at(iPath);
            auto nSegments = path.size();
            std::vector<std::pair<double, double>> polyPath;
            polyPath.reserve(nSegments + 1);
            size_t i = 0;
            for (const auto &pi : path)
            {
                auto point = pi.getStartPoint();
                auto xi = point.getPositionInX();
                auto zi = point.getPositionInZ();
                polyPath.push_back(std::pair {xi, zi});
                if (i == nSegments - 1)
                {
                    point = pi.getEndPoint();
                    xi = point.getPositionInX();
                    zi = point.getPositionInZ();
                    polyPath.push_back(std::pair {xi, zi}); 
                }
                i = i + 1;
            }
            lines.push_back(std::move(polyPath));
        }
        vtkWriter.write(lines);
    }
    vtkWriter.close();
}

/// Get the ray paths
std::vector<Path2D> GradientTracer2D::getRayPaths() const
{
    if (!haveRayPaths()){throw std::runtime_error("Ray paths not computed");}
    return pImpl->mPaths;
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<double> &solver);
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<float> &solver);
