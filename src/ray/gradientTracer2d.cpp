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
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/abstractBaseClass/solver2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/station2d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX::Ray;

namespace
{

enum class Wall : int
{
    Top = 0,
    Right = 1,
    Bottom = 2,
    Left = 3,
    Unknown =-1
};

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

/*
template<typename T>
T bilinear(const T x, const T z,
           const T x1, const T x2, 
           const T z1, const T z2, 
           const T f00, const T f01,
           const T f10, const T f11)
{
    T dxi = static_cast<T> (1./(x2 - x1));
    T dzi = static_cast<T> (1./(z2 - z1));
    return ::bilinear(x, z, 
                      x1, x2,
                      z1, z2,
                      f00, f01,
                      f10, f11,
                      dxi, dzi);
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
/// Note, the four (left, right, bottom, top) cell faces are:
///   F_l = x_0
///   F_r = x_1
///   F_b = z_0 
///   F_t = z_1
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
/*
template<typename T>
void intersect(const int iCurrentCellX, const int iCurrentCellZ,
               const int nCellX, const int nCellZ,
               const T ox, const T oz,
               const T gx, const T gz,
               const T dx, const T dz,
               T *ex, T *ez)
{
    const double tol{1.e-2}; // 2 cm
    *ex = ox;
    *ez = oz;
    T ti{1.e10};
    T tj{1.e10};
    T gxUse = gx;
    T gzUse = gz;
    // The easy case is when we are in the cell.
    int ix0 = iCurrentCellX;
    int ix1 = ix0 + 1;
    int iz0 = iCurrentCellZ;
    int iz1 = iz0 + 1;
    // It's more complicated when we are on a boundary:
    // Left cell boundry
    if (std::abs(ox - iCurrentCellX*dx) < tol)
    {
        // Make sure gradient doesn't take us out of bounds
        if (iCurrentCellX == 0)
        {
            if (gx < 0){gxUse = 0;}
            ix0 = 0;
            ix1 = 1;
        }
        else
        {
            if (gx < 0)
            {
                ix0 = iCurrentCellX - 1;
                ix1 = iCurrentCellX;
            }
            else
            {
                ix0 = iCurrentCellX;
                ix1 = iCurrentCellX + 1;
            }
        }
    }
    // Right cell boundary
    else if (std::abs(ox - (iCurrentCellX + 1)*dx) < tol)
    {
        if (iCurrentCellX == nCellX - 1)
        {
            if (gx > 0){gxUse = 0;}
            ix0 = nCellX - 1;
            ix1 = nCellX;
        }
        else
        {
            if (gx > 0)
            {
                ix0 = iCurrentCellX + 1;
                ix1 = iCurrentCellX + 2;
            }
            else
            {
                ix0 = iCurrentCellX;
                ix1 = iCurrentCellX + 1;
            }
        }
    }
    if (std::abs(oz - iCurrentCellZ*dz) < tol)
    {
        if (iCurrentCellZ == 0)
        {
            if (gz < 0){gzUse = 0;}
            iz0 = 0;
            iz1 = 1;
        }
        else
        {
            if (gz < 0)
            {
                iz0 = iCurrentCellZ - 1;
                iz1 = iCurrentCellZ;
            }
            else
            {
                iz0 = iCurrentCellZ;
                iz1 = iCurrentCellZ + 1;
            }
        }
    }
    else if (std::abs(oz - (iCurrentCellZ + 1)*dz) < tol)
    {
        if (iCurrentCellZ == nCellZ - 1)
        {
            if (gz > 0){gzUse = 0;}
            iz0 = nCellZ - 1;
            iz1 = nCellZ;    
        }
        else
        {
            if (gz > 0)
            {
                iz0 = iCurrentCellZ + 1;
                iz1 = iCurrentCellZ + 2;
            }
            else
            {
                iz0 = iCurrentCellZ;
                iz1 = iCurrentCellZ + 1;
            }
        }
    }
    // Normalize vectors
    double gNorm = std::hypot(gxUse, gzUse);
    gxUse = gxUse/gNorm;
    gzUse = gzUse/gNorm;
    // Figure out the horizontal then vertical planes
    T z0 = oz + iz0*dz;
    T z1 = oz + iz1*dz;
    T x0 = ox + ix0*dx;
    T x1 = ox + ix0*dx;
std::cout << "gUse: " << gxUse << " " << gzUse << std::endl;
std::cout << ix0 << " " << ix1 << " " << iz0 << " " << iz1 << std::endl;
    /// Existence
    if (gxUse > 0)
    {
        ti = (x1 - ox)/gxUse;
    }
    else if (gxUse < 0)
    {
        ti = (x0 - ox)/gxUse;
    }
    if (gzUse > 0)
    {
        tj = (z1 - oz)/gzUse;
    }
    else if (gz > 0)
    {
        tj = (z0 - oz)/gzUse;
    }
    // Make my lines
    T fxi = ox + gxUse*ti;
    T fzi = oz + gzUse*ti;
    T fxj = ox + gxUse*tj;
    T fzj = oz + gzUse*tj;
std::cout << "next point i: " << fxi << " " << fzi << std::endl;
std::cout << "next point j: " << fxj << " " << fzj << std::endl;
    // Distance squared is good enough
    T li = (fxi - ox)*(fxi - ox) + (fzi - oz)*(fzi - oz);
    T lj = (fxj - ox)*(fxj - ox) + (fzj - oz)*(fzj - oz);
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

template<typename T>
void computeNextPointInPath(const int iCurrentCellX, const int iCurrentCellZ,
                            const double xRay0, const double zRay0,
                            const double dx, const double dz,
                            const T *gradientPtr,
                            const int nGridX, const int nGridZ,
                            int *iNextCellX, int *iNextCellZ,
                            double *xRay1, double *zRay1)
{
    *iNextCellX = iCurrentCellX;
    *iNextCellZ = iCurrentCellZ;
    auto i00 = 2*::gridToIndex(nGridX, iCurrentCellX,     iCurrentCellZ);
    auto i10 = 2*::gridToIndex(nGridX, iCurrentCellX + 1, iCurrentCellZ);
    auto i11 = 2*::gridToIndex(nGridX, iCurrentCellX + 1, iCurrentCellZ + 1); 
    auto i01 = 2*::gridToIndex(nGridX, iCurrentCellX,     iCurrentCellZ + 1); 
    auto gx00 = static_cast<double> (gradientPtr[i00]);
    auto gz00 = static_cast<double> (gradientPtr[i00 + 1]);
    auto gx10 = static_cast<double> (gradientPtr[i10]);
    auto gz10 = static_cast<double> (gradientPtr[i10 + 1]);
    auto gx01 = static_cast<double> (gradientPtr[i01]);
    auto gz01 = static_cast<double> (gradientPtr[i01 + 1]);
    auto gx11 = static_cast<double> (gradientPtr[i11]);
    auto gz11 = static_cast<double> (gradientPtr[i11 + 1]); 
    // Interpolate gradient.  Note, we reduce to unit cell [0, dx] x [0, dz].
    double x0 = iCurrentCellX*dx;
    double x1 = (iCurrentCellX + 1)*dx;
    double z0 = iCurrentCellZ*dz;
    double z1 = (iCurrentCellZ + 1)*dz;
    constexpr double zero{0};
    auto gx = ::bilinear(xRay0 - x0, zRay0 - z0,
                         zero, dx, zero, dz,
                         gx00, gx01, gx10, gx11,
                         1./dx, 1./dz);
    auto gz = ::bilinear(xRay0 - x0, zRay0 - z0,
                         zero, dx, zero, dz,
                         gz00, gz01, gz10, gz11,
                         1./dx, 1./dz);
    // March down the gradient
    gx =-gx;
    gz =-gz;
    // Figure out where next to advance ray
std::cout << xRay0 << " " << zRay0 << std::endl;
std::cout << gx << " " << gz << std::endl;
    ::intersect(iCurrentCellX, iCurrentCellZ,
                nGridX - 1, nGridZ - 1,
                xRay0, zRay0, gx, gz, dx, dz, xRay1, zRay1);
}
*/

}

class GradientTracer2D::GradientTracer2DImpl
{
public:
    EikonalXX::Geometry2D mGeometry;
    std::vector<EikonalXX::Station2D> mStations;
    std::vector<EikonalXX::Ray::Path2D> mPaths;
    bool mInitialized{false};
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
}

/// Have stations?
bool GradientTracer2D::haveStations() const noexcept
{
    return !pImpl->mStations.empty();
}

/// Set the geometry
void GradientTracer2D::initialize(const EikonalXX::Geometry2D &geometry)
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
    pImpl->mGeometry = geometry;
    pImpl->mInitialized = true;
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
    // Grid properties
    int nGridX = pImpl->mGeometry.getNumberOfGridPointsInX();
    int nGridZ = pImpl->mGeometry.getNumberOfGridPointsInZ();
    double dx = pImpl->mGeometry.getGridSpacingInX();
    double dz = pImpl->mGeometry.getGridSpacingInZ();
    double dxi = 1./dx;
    double dzi = 1./dz;
    double x0 = pImpl->mGeometry.getOriginInX();
    double z0 = pImpl->mGeometry.getOriginInZ();
    double stepLength = 0.1*std::min(dx, dz);
    // Source properties
    auto source = solver.getSource();
    auto iSourceCellX = source.getCellInX();
    auto iSourceCellZ = source.getCellInZ();
    auto iSourceCell = source.getCell();
    auto xs = source.getOffsetInX();
    auto zs = source.getOffsetInZ();
std::cout << source << std::endl;
    // Set the source point (we'll always use this when tracing)
    Point2D sourcePoint{x0 + xs, z0 + zs}; 
    const T* gradientPtr = solver.getTravelTimeGradientFieldPointer();
    // Loop on receivers
    for (const auto &station : pImpl->mStations)
    {
        auto iStationCellX = station.getCellInX();
        auto iStationCellZ = station.getCellInZ();
        // Define the station piont (we'll also always use this when tracing)
        auto xr = station.getOffsetInX();
        auto zr = station.getOffsetInZ();
        Point2D stationPoint {x0 + xr, z0 + zr};
        auto iStationCell = station.getCell();
        // Deal with an algorithm breakdown where the source/station are in
        // the same cell
        if (iSourceCell == iStationCell)
        {
            Segment2D segment;
            segment.setStartAndEndPoint(std::pair {sourcePoint, stationPoint});
            segment.setSlowness(
                solver.getSlowness(station.getCellInX(), station.getCellInZ()));
            segment.setVelocityModelCellIndex(iStationCell); 
            continue; 
        }
        // Business as usual - march down the gradient
        double xRay0 = xr;
        double zRay0 = zr;
        double xRay1 = 0;
        double zRay1 = 0;
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
        bool mConverged{false};
        const int maxSegments{2500};
        for (int iSegment = 0; iSegment < maxSegments; ++iSegment)
        {
            auto iCellX = static_cast<int> (xRay0/dx);
            auto iCellZ = static_cast<int> (zRay0/dz);
            auto gx00 = gx000;
            auto gz00 = gz000;
            auto gx10 = gx100;
            auto gz10 = gz100;
            auto gx01 = gx010;
            auto gz01 = gz010;
            auto gx11 = gx110;
            auto gz11 = gz110; 
            if (iCellX != iCellX0 || iCellZ != iCellZ0)
            {
                auto i00 = 2*::gridToIndex(nGridX, iCellX,     iCellZ);
                auto i10 = 2*::gridToIndex(nGridX, iCellX + 1, iCellZ);
                auto i11 = 2*::gridToIndex(nGridX, iCellX + 1, iCellZ + 1);
                auto i01 = 2*::gridToIndex(nGridX, iCellX,     iCellZ + 1);
                gx00 = static_cast<double> (gradientPtr[i00]);
                gz00 = static_cast<double> (gradientPtr[i00 + 1]);
                gx10 = static_cast<double> (gradientPtr[i10]);
                gz10 = static_cast<double> (gradientPtr[i10 + 1]);
                gx01 = static_cast<double> (gradientPtr[i01]);
                gz01 = static_cast<double> (gradientPtr[i01 + 1]);
                gx11 = static_cast<double> (gradientPtr[i11]);
                gz11 = static_cast<double> (gradientPtr[i11 + 1]); 
            }
            auto gx = ::bilinear(xRay0, zRay0,
                                 iCellX*dx, (iCellX + 1)*dx,
                                 iCellZ*dz, (iCellZ + 1)*dz,
                                 gx00, gx01,
                                 gx10, gx11,
                                 dxi, dzi);
            auto gz = ::bilinear(xRay0, zRay0,
                                 iCellX*dx, (iCellX + 1)*dx,
                                 iCellZ*dz, (iCellZ + 1)*dz,
                                 gz00, gz01,
                                 gz10, gz11,
                                 dxi, dzi);
            auto gNorm = std::hypot(gx, gz);
            auto gStep = stepLength/gNorm;
            // March a small step opposite direction of gradient
            xRay1 = xRay0 - gx*gStep;
            zRay1 = zRay0 - gz*gStep;
            std::cout << xRay1 << " " << zRay1 << std::endl;
            // Did we converge?
            if (std::abs(iCellX - iSourceCellX) < 2 &&
                std::abs(iCellZ - iSourceCellZ) < 2)
            {
                mConverged = true;
                break;
            }
            // Update
            xRay0 = xRay1;
            zRay0 = zRay1;
            iCellX0 = iCellX;
            iCellZ0 = iCellZ;
            gx000 = gx00;
            gz000 = gz00;
            gx100 = gx10;
            gz100 = gz10;
            gx010 = gx01;
            gz010 = gz01;
            gx110 = gx11;
            gz110 = gz11;
/*
            //  
            computeNextPointInPath(iCurrentCellX, iCurrentCellZ,
                                   xRay0, zRay0,
                                   dx, dz,
                                   gradientPtr,
                                   nGridX, nGridZ,
                                   &iNextCellX, &iNextCellZ,
                                   &xRay1, &zRay1);
            Point2D p0{xRay0 + x0, zRay0 + z0};
            Point2D p1{xRay1 + x0, zRay1 + z0};
std::cout << xRay0 + x0 << " " << zRay0 + z0 << std::endl;
std::cout << xRay1 + x0 << " " << zRay1 + z0 << std::endl;
            Segment2D segment;
            segment.setStartAndEndPoint(std::pair {p0, p1});
            // Are we at the source?
getchar();
*/
        }
        if (mConverged){std::cout << "Converged!" << std::endl;}
getchar();
        if (!mConverged)
        {
            std::cerr << "Ray for station: " << station.getName()
                      << " did not converge to source" << std::endl;
        }
    }
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<double> &solver);
template void EikonalXX::Ray::GradientTracer2D::trace(
    const EikonalXX::AbstractBaseClass::ISolver2D<float> &solver);
