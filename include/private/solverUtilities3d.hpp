#ifndef PRIVATE_SOLVERUTILITIES3D_HPP
#define PRIVATE_SOLVERUTILITIES3D_HPP
#include <CL/sycl.hpp>
#include <cmath>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#include <tbb/tbb.h>
#include <tbb/cache_aligned_allocator.h>
#include "eikonalxx/enums.hpp"
#include "eikonalxx/graph3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "private/grid.hpp"
#include "private/pad.hpp"

#define BOUNDARY_NODE -1
#define SOURCE_NODE 0
#define UPDATE_NODE 1
namespace
{

template<EikonalXX::SweepNumber3D E>
void getSweepFiniteDifferenceSigns(int *ixShift, int *iyShift, int *izShift,
                                   int *signX, int *signY, int *signZ);
template<EikonalXX::SweepNumber3D E>
void gridToSurroundingSlowness(
    const int ix, const int iy, const int iz,
    const int nCellX, const int nCellY, const int nCellZ,
    int *iCell0, int *iCell1, int *iCell2, int *iCell3,
    int *iCell4, int *iCell5, int *iCell7);

using namespace EikonalXX;

/// @brief Container for the slownesses in a level for the 2D solver.
template<class T>
struct SweepSlowness3D
{
    /// Allocates space to hold slownesses for the n nodes in the sweep.
    void allocate(int n)
    {
        clear();
        if (n <= 0){throw std::invalid_argument("n must be positive");}
        nNodes = n; 
        s0.resize(static_cast<size_t> (nNodes), 0);
        s1.resize(static_cast<size_t> (nNodes), 0);
        s2.resize(static_cast<size_t> (nNodes), 0);
        s3.resize(static_cast<size_t> (nNodes), 0);
        s4.resize(static_cast<size_t> (nNodes), 0);
        s5.resize(static_cast<size_t> (nNodes), 0);
        s7.resize(static_cast<size_t> (nNodes), 0); 
    }    
    /// This sets all the slownesses in the level to 0
    void zero() noexcept
    {
        if (nNodes > 0) 
        {
            std::fill(s0.begin(), s0.end(), 0);
            std::fill(s1.begin(), s1.end(), 0);
            std::fill(s2.begin(), s2.end(), 0);
            std::fill(s3.begin(), s3.end(), 0);
            std::fill(s4.begin(), s4.end(), 0);
            std::fill(s5.begin(), s5.end(), 0);
            std::fill(s7.begin(), s7.end(), 0);
        }    
    }    
    /// Gets the min slowness in the sweep
    T getMinimumValue() const noexcept    
    {    
        T sMin = 0;
        if (nNodes > 0)
        {
            sMin = *std::min_element(s0.begin(), s0.end());
            sMin = std::min(sMin, *std::min_element(s1.begin(), s1.end()));
            sMin = std::min(sMin, *std::min_element(s2.begin(), s2.end()));
            sMin = std::min(sMin, *std::min_element(s3.begin(), s3.end()));
            sMin = std::min(sMin, *std::min_element(s4.begin(), s4.end()));
            sMin = std::min(sMin, *std::min_element(s5.begin(), s5.end()));
            sMin = std::min(sMin, *std::min_element(s7.begin(), s7.end()));
        }
        return sMin;
    }    
    /// Releases memory.
    void clear() noexcept
    {
        s0.clear();
        s1.clear();
        s2.clear();
        s3.clear();
        s4.clear();
        s5.clear();
        s7.clear();
        nNodes = 0; 
    }    
    /// The slowness (s/m) in the home cell.
    /// This is an array whose dimension is at least [nNodes].
    std::vector<T, tbb::cache_aligned_allocator<T>> s0;
    /// The slowness (s/m) in the cell to the right or left of the home cell.
    /// This is an array whose dimension is at least [nNodes].
    std::vector<T, tbb::cache_aligned_allocator<T>> s1;
    /// The slowness (s/m) in the cell diagonal but in the same plane as the
    /// home cell.
    std::vector<T, tbb::cache_aligned_allocator<T>> s2;
    /// The slowness (s/m) in the cell in front of or in back of the home cell.
    /// This is an array whose dimension is at least [nNodes].
    std::vector<T, tbb::cache_aligned_allocator<T>> s3;
    /// The slowness (s/m) in the cell above or below the home cell.
    /// This is an array whose dimension is at least [nNodes].
    std::vector<T, tbb::cache_aligned_allocator<T>> s4;
    /// The slowness (s/m) in the cell above or below cell s1.
    std::vector<T, tbb::cache_aligned_allocator<T>> s5;
    /// The slowness (s/m) in the cell above or below cell s3.
    std::vector<T, tbb::cache_aligned_allocator<T>> s7;
    /// The number of nodes in the sweep.
    int nNodes = 0; 
};

///--------------------------------------------------------------------------///
///                             Analytical Solutions                         ///
///--------------------------------------------------------------------------///
/// @brief Computes the travel time from a source to a grid point in a 
///        homogeneous media.
/// @param[in] ix             The x grid index.
/// @param[in] iy             The y grid index.
/// @param[in] iz             The z grid index.
/// @param[in] dx             The grid spacing in x in meters.
/// @param[in] dy             The grid spacing in y in meters.
/// @param[in] dz             The grid spacing in z in meters.
/// @param[in] xSourceOffset  The source's x offset from the origin in meters.
///                           This is the source x location - model x origin.
/// @param[in] ySourceOffset  The source's y offset from the origin in meters.
///                           This is the source y location - model y origin.
/// @param[in] zSourceOffset  The source's z offset from the origin in meters.
///                           This is the source z location - model z origin.
/// @param[in] slowness       The slowness at the source in s/m. 
/// @result The travel time from the source to the grid point in seconds.
#pragma omp declare simd uniform(dx, dy, dz, xSourceOffset, ySourceOffset, zSourceOffset)
template<class T>
[[nodiscard]]
T computeAnalyticalTravelTime(const int ix, const int iy, const int iz,
                              const T dx, const T dy, const T dz,
                              const T xSourceOffset,
                              const T ySourceOffset,
                              const T zSourceOffset,
                              const T slowness)
{
    auto x = ix*dx;
    auto y = iy*dy;
    auto z = iz*dz;
    auto deltaX = x - xSourceOffset;
    auto deltaY = y - ySourceOffset;
    auto deltaZ = z - zSourceOffset;
    auto distance = sycl::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return slowness*distance;
}
/// @brief Computes the travel time from a source to a grid point at as
///        well as the derivatives at the grid point in a homogeneous media.
/// @param[in] ix             The x grid index.
/// @param[in] iy             The y grid index.
/// @param[in] iz             The z grid index.
/// @param[in] dx             The grid spacing in x in meters.
/// @param[in] dy             The grid spacing in y in meters.
/// @param[in] dz             The grid spacing in z in meters.
/// @param[in] xSourceOffset  The source's x offset from the origin in meters.
///                           This is the source x location - model x origin.
/// @param[in] ySourceOffset  The source's y offset from the origin in meters.
///                           This is the source y location - model y origin.
/// @param[in] zSourceOffset  The source's z offset from the origin in meters.
///                           This is the source z location - model z origin.
/// @param[in] slowness       The slowness at the source in s/m.
/// @param[out] t             The travel time from the source to the grid
///                           point in seconds.
/// @param[out] dtdx          The derivative of the travel time w.r.t. x
///                           at the grid point in s/m. 
/// @param[out] dtdy          The derivative of the travel time w.r.t. y
///                           at the grid point in s/m.
/// @param[out] dtdz          The derivative of the travel time w.r.t. z
///                           at the grid point in s/m.
#pragma omp declare simd uniform(dx, dy, dz, xSourceOffset, ySourceOffset, zSourceOffset)
template<class T>
void computeAnalyticalTravelTime(const int ix, const int iy, const int iz,
                                 const T dx, const T dy, const T dz,
                                 const T xSourceOffset,
                                 const T ySourceOffset,
                                 const T zSourceOffset,
                                 const T slowness,
                                 T *t, T *dtdx, T *dtdy, T *dtdz)
{
    auto x = ix*dx;
    auto y = iy*dy;
    auto z = iz*dz;
    auto deltaX = x - xSourceOffset;
    auto deltaY = y - ySourceOffset;
    auto deltaZ = z - zSourceOffset;
    auto distance = sycl::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    *t = slowness*distance;
    // limit approaching from right of x/sqrt(x) is 0
    *dtdx = 0; //std::numeric_limits<T>::quiet_NaN();
    *dtdy = 0; //std::numeric_limits<T>::quiet_NaN();
    *dtdz = 0; //std::numeric_limits<T>::quiet_NaN();
    if (distance > 0) 
    {    
        *dtdx = (slowness*deltaX)/distance;
        *dtdy = (slowness*deltaY)/distance;
        *dtdz = (slowness*deltaZ)/distance; 
    }    
}
///--------------------------------------------------------------------------///
///                              Finite Differences                          ///
///--------------------------------------------------------------------------///

template<typename T>
[[nodiscard]]
T finiteDifference(const int sphericalRadius,
                   const T huge,
                   const T h,
                   const int ix, const int iy, const int iz,
                   const int signX, const int signY, const int signZ,
                   const int ixShift, const int iyShift, const int izShift,
                   const int iSrcX, const int iSrcY, const int iSrcZ,
                   const T xSourceOffset,
                   const T ySourceOffset,
                   const T zSourceOffset,
                   const T s0, const T s1, const T s2, const T s3,
                   const T s4, const T s5, const T s7,
                   const T t1, const T t2, const T t3,
                   const T t4, const T t5, const T t6, const T t7)
{
    // Default result
    T tUpd = huge;
    // Some mins I'll need for the refracted plane operators
    auto minS0S3 = sycl::fmin(s0, s3); 
    auto minS0S1 = sycl::fmin(s0, s1);
    auto minS0S4 = sycl::fmin(s0, s4);
    // 1D operators critically refracted operators
    auto t1d1 = t1 + h*sycl::fmin(minS0S3, sycl::fmin(s4, s7)); 
    auto t1d2 = t3 + h*sycl::fmin(minS0S4, sycl::fmin(s1, s5));
    auto t1d3 = t4 + h*sycl::fmin(minS0S1, sycl::fmin(s2, s3));
    auto t1d = sycl::min(sycl::fmin(t1d1, t1d2), t1d3);
    // Cartesian
    T t2d = 0;
    if (sycl::abs(iSrcX - ix) > sphericalRadius ||
        sycl::abs(iSrcY - iy) > sphericalRadius ||
        sycl::abs(iSrcZ - iz) > sphericalRadius)
    {
        // 2D critically refracted operators
        T dTxXZ =-t1 + t4 - t5;
        T dTzXZ =-t4 + t1 - t5;
        T dtCrossXZ = t1 - t4;
        T hs03 = h*minS0S3;
        T detXZ = 2*(hs03*hs03) - dtCrossXZ*dtCrossXZ;
        T t2d1 = huge;
        if (detXZ > 0 && t4 <= t1 + hs03 && t1 <= t4 + hs03)
        {
            t2d1 = t5 + sycl::sqrt(detXZ);
        }
        t2d = sycl::fmin(t1d, t2d1);

        T dTyYZ =-t3 + t4 - t7;
        T dTzYZ =-t4 + t3 - t7;
        T dtCrossYZ = t3 - t4;
        T hs01 = h*minS0S1;
        T detYZ = 2*(hs01*hs01) - dtCrossYZ*dtCrossYZ;
        T t2d2 = huge;
        if (detYZ > 0 && t3 <= t4 + hs01 && t4 <= t3 + hs01)
        {
            t2d2 = t7 + sycl::sqrt(detYZ);
        }
        t2d = sycl::fmin(t2d, t2d2);

        T dTxXY =-t1 + t3 - t2;
        T dTyXY =-t3 + t1 - t2;
        T dtCrossXY = t1 - t3;
        T hs04 = h*minS0S4;
        T detXY = 2*(hs04*hs04) - dtCrossXY*dtCrossXY;
        T t2d3 = huge;
        if (detXY > 0 && t3 <= t1 + hs04 && t1 <= t3 + hs04)
        {
            t2d3 = t2 + sycl::sqrt(detXY);
        }
        t2d = sycl::fmin(t2d, t2d3);

        // 3D operator
        T t3d = huge;

    }
    else
    {
        T t0, dtdx, dtdy, dtdz;
        computeAnalyticalTravelTime(ix, iy, iz,
                                    h, h, h,
                                    xSourceOffset, ySourceOffset, zSourceOffset,
                                    s0,
                                    &t0, &dtdx, &dtdy, &dtdz);
        T tau1 = t1 - computeAnalyticalTravelTime(ix + ixShift,
                                                  iy,
                                                  iz, 
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);  
        T tau2 = t2 - computeAnalyticalTravelTime(ix + ixShift,
                                                  iy + iyShift,
                                                  iz,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);
        T tau3 = t3 - computeAnalyticalTravelTime(ix,
                                                  iy + iyShift,
                                                  iz,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);
        T tau4 = t4 - computeAnalyticalTravelTime(ix,
                                                  iy, 
                                                  iz + izShift,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);
        T tau5 = t5 - computeAnalyticalTravelTime(ix + ixShift,
                                                  iy,
                                                  iz + izShift,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);  
        T tau6 = t6 - computeAnalyticalTravelTime(ix + ixShift,
                                                  iy + iyShift,
                                                  iz + izShift,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);
        T tau7 = t7 - computeAnalyticalTravelTime(ix,
                                                  iy + iyShift,
                                                  iz + izShift,
                                                  h, h, h,
                                                  xSourceOffset,
                                                  ySourceOffset,
                                                  zSourceOffset,
                                                  s0);
    }
    return tUpd;
}


///--------------------------------------------------------------------------///
///                              Grid Calculations                           ///
///--------------------------------------------------------------------------///

/// @brief Gets the sweep finite difference signs and shifts.
/// @param[out] ixShift  The shift in the x grid point when computing the
///                      analytic travel times and derivatives.
/// @param[out] iyShift  The shift in the y grid point when computing the
///                      analytic travel times and derivatives.
/// @param[out] izShift  The shift in the z grid point when computing the
///                      analytic travel times and derivatives.
/// @param[out] signX    The sign on the derivative in x.  This makes the
///                      derivative sign consistent with the sweep direction.
/// @param[out] signY    The sign on the derivative in y.  This makes the
///                      derivative sign consistent with the sweep direction.
/// @param[out] signZ    The sign on the derivative in z.  This makes the
///                      derivative sign consistent with the sweep direction.
template<EikonalXX::SweepNumber3D E>
void getSweepFiniteDifferenceSigns(int *ixShift, int *iyShift, int *izShift,
                                   int *signX, int *signY, int *signZ)
{
    if constexpr(E == EikonalXX::SweepNumber3D::SWEEP1)
    {
        *ixShift =-1;
        *iyShift =-1;
        *izShift =-1;
        *signX = 1;
        *signY = 1;
        *signZ = 1; 
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP2)
    {
        *ixShift =+1;
        *iyShift =-1;
        *izShift =-1;
        *signX =-1;
        *signY = 1;
        *signZ = 1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP3)
    {
        *ixShift =-1;
        *iyShift =+1;
        *izShift =-1;
        *signX = 1;
        *signY =-1;
        *signZ = 1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP4)
    {
        *ixShift =+1;
        *iyShift =+1;
        *izShift =-1;
        *signX =-1;
        *signY =-1;
        *signZ = 1; 
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP5) // Sweep 1; flip z
    {
        *ixShift =-1;
        *iyShift =-1;
        *izShift =+1;
        *signX = 1;
        *signY = 1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP6) // Sweep 2; flip z
    {
        *ixShift =+1;
        *iyShift =-1;
        *izShift =+1;
        *signX =-1;
        *signY = 1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP7) // Sweep 3; flip z
    {
        *ixShift =-1;
        *iyShift =+1;
        *izShift =+1;
        *signX = 1;
        *signY =-1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::SWEEP8) // Sweep 4 ; flip z
    {
        *ixShift =+1;
        *iyShift =+1;
        *izShift =+1;
        *signX =-1;
        *signY =-1;
        *signZ =-1;
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

/// @brief Gets the loop limits for the fast sweeping method i.e.,:
///        for (int iz = iz0; iz != iz1; iz = iz + izDir)
///            for (int iy = iy0; iy != iy1; iy = iy + iyDir)
///                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
/// @param[in] sweep   The sweep number.
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridY  The number of grid points in y.
/// @param[in] nGridZ  The number of grid points in z.
/// @param[out] ix0    The first x grid index in the loop.
/// @param[out] ix1    The stopping x grid index (exclusive).
/// @param[out] iy0    The first y grid index in the loop.
/// @param[out] iy1    The stopping y grid index (exclusive).
/// @param[out] iz0    The first z grid index in the loop.
/// @param[out] iz1    The stopping z grid index (exclusive).
/// @param[out] ixDir  The x loop variable update increment (+1 or -1).
/// @param[out] iyDir  The y loop variable update increment (+1 or -1).
/// @param[out] izDir  The z loop variable update increment (+1 or -1).
template<EikonalXX::SweepNumber3D E>
void getLoopLimits(const int nGridX, const int nGridY, const int nGridZ,
                   int *ix0, int *iy0, int *iz0,
                   int *ix1, int *iy1, int *iz1,
                   int *ixDir, int *iyDir, int *izDir)
{
    if constexpr (E == SweepNumber3D::SWEEP1)
    {
        *ix0 = 1;
        *ix1 = nGridX;
        *iy0 = 1;
        *iy1 = nGridY;
        *iz0 = 1;
        *iz1 = nGridZ;
        *ixDir = 1;
        *iyDir = 1;
        *izDir = 1;
    }    
    else if constexpr (E == SweepNumber3D::SWEEP2)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iy0 = 1;
        *iy1 = nGridY;
        *iz0 = 1;
        *iz1 = nGridZ;
        *ixDir =-1;
        *iyDir = 1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP3)
    {
        *ix0 = 1; 
        *ix1 = nGridX;
        *iy0 = nGridY - 2; 
        *iy1 =-1; 
        *iz0 = 1;  
        *iz1 = nGridZ;
        *ixDir = 1; 
        *iyDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP4)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iy0 = nGridY - 2;
        *iy1 =-1;
        *iz0 = 1;
        *iz1 = nGridZ;
        *ixDir =-1;
        *iyDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP5)
    {
        *ix0 = 1;
        *ix1 = nGridX;
        *iy0 = 1;
        *iy1 = nGridY;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir = 1;
        *iyDir = 1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP6)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iy0 = 1;
        *iy1 = nGridY;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir =-1;
        *iyDir = 1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP7)
    {
        *ix0 = 1;
        *ix1 = nGridX;
        *iy0 = nGridY - 2;
        *iy1 =-1;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir = 1;
        *iyDir =-1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::SWEEP8)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iy0 = nGridY - 2;
        *iy1 =-1;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir =-1;
        *iyDir =-1;
        *izDir =-1;
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

/// @brief Given an (ix, iy, iz) grid point this permutes to the corresponding
///        (jx, jy, jz) in the sweep.
/// @param[in] ix   The grid index in x.
/// @param[in] iy   The grid index in y.
/// @param[in] iz   The grid index in z.
/// @param[in] nx   The number of grid points in x.
/// @param[in] ny   The number of grid points in y.
/// @param[in] nz   The number of grid points in z.
/// @param[out] jx  The corresponding x grid point in the sweep.
/// @param[out] jy  The corresponding y grid point in the sweep.
/// @param[out] jz  The corresponding z grid point in the sweep.
#pragma omp declare simd uniform(nx, ny, nz)
template<EikonalXX::SweepNumber3D E>
void permuteGrid(const int ix, const int iy, const int iz,
                 const int nx, const int ny, const int nz,
                 int *jx, int *jy, int *jz)
{
    if constexpr (E == SweepNumber3D::SWEEP1)
    {
        *jx = ix;
        *jy = iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP2)
    {
        *jx = nx - 1 - ix;
        *jy = iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP3)
    {
        *jx = ix;
        *jy = ny - 1 - iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP4)
    {
        *jx = nx - 1 - ix;
        *jy = ny - 1 - iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP5)
    {
        *jx = ix;
        *jy = iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP6)
    {
        *jx = nx - 1 - ix;
        *jy = iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP7)
    {
        *jx = ix;
        *jy = ny - 1 - iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::SWEEP8)
    {
        *jx = nx - 1 - ix;
        *jy = ny - 1 - iy;
        *jz = nz - 1 - iz;
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

/// @brief Converts a grid point (ix, iy, iz) in a sweep to the corresponding
///        level and index.
/// @param[in] ix      The grid index in x.
/// @param[in] iy      The grid index in y.
/// @param[in] iz      The grid index in z.
/// @param[in] nx      The number of grid points in x.
/// @param[in] ny      The number of grid points in y.
/// @param[in] nz      The number of grid points in z.
/// @param[out] level  The corresponding level, which is a plane defined by
///                    constant = ix + iy + iz.
/// @param[out] indx   The corresponding indx in the level.  Realize, this does
///                    not necessarily begin at 0.  To obtain this information
///                    use \c getLevelStartStopIndices() and subtract by i0.
/// @note This is the inverse operation of \c sweepLevelIndexToGrid().
#pragma omp declare simd uniform(nx, ny, nz)
template<EikonalXX::SweepNumber3D E>
void gridToLevel(const int ix, const int iy, const int iz,
                 const int nx, const int ny, const int nz,
                 int *level)
{
    int jx, jy, jz;
    permuteGrid<E>(ix, iy, iz, nx, ny, nz, &jx, &jy, &jz);
    *level = jx + jy + jz;
}

/// @param[out] iCell0  The index of the slowness field corresponding to
///                     the home cell.
/// @param[out] iCell1  The index of the slowness field corresponding to
///                     the cell to the left or right of the home cell.
/// @param[out] iCell2  The index of the slowness field corresponding to
///                     the cell adjacent from but in the same plane as
///                     the home cell.
/// @param[out] icell3  The index of the slowness field corresponding to
///                     the cell in front of or in back of the home cell.
/// @param[out] iCell4  The index of the slowness field corresponding to
///                     the cell above or below the home cell.
/// @param[out] iCell5  The index of the slowness field corresponding to
///                     the cell to the lower right or upper right of
///                     the home cell.
/// @param[out] iCell7  The index of the slowness field corresponding to
///                     the cell to the 
#pragma omp declare simd uniform(nCellX, nCellY, nCellZ)
template<EikonalXX::SweepNumber3D E>
void gridToSurroundingSlownessIndices(
    const int ix, const int iy, const int iz,
    const int nCellX, const int nCellY, const int nCellZ,
    int *iCell0, int *iCell1, int *iCell2, int *iCell3,
    int *iCell4, int *iCell5, int *iCell7)
{
    if constexpr (E == SweepNumber3D::SWEEP1)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix - 1));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate + 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy - 1));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate + 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz - 1));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate + 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate); 
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP2)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate - 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy - 1));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate + 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz - 1));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate + 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP3)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix - 1));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate + 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate - 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz - 1));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate + 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP4)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate - 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate - 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz - 1));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate + 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP5)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix - 1));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate + 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy - 1));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate + 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate - 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP6)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate - 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy - 1));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate + 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate - 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP7)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix - 1));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate + 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate - 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate - 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
    else if constexpr (E == SweepNumber3D::SWEEP8)
    {
        int ixUpdate = std::max(0, std::min(nCellX - 1, ix));
        int ixOffset = std::max(0, std::min(nCellX - 1, ixUpdate - 1));
        int iyUpdate = std::max(0, std::min(nCellY - 1, iy));
        int iyOffset = std::max(0, std::min(nCellY - 1, iyUpdate - 1));
        int izUpdate = std::max(0, std::min(nCellZ - 1, iz));
        int izOffset = std::max(0, std::min(nCellZ - 1, izUpdate - 1));
        *iCell0 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izUpdate);
        *iCell1 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izUpdate);
        *iCell2 = gridToIndex(nCellX, nCellY, ixOffset, iyOffset, izUpdate);
        *iCell3 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izUpdate);
        *iCell4 = gridToIndex(nCellX, nCellY, ixUpdate, iyUpdate, izOffset);
        *iCell5 = gridToIndex(nCellX, nCellY, ixOffset, iyUpdate, izOffset);
        *iCell7 = gridToIndex(nCellX, nCellY, ixUpdate, iyOffset, izOffset);
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

/// @brief Converts a grid index to the surrounding travel time indices 
///        used in the finite difference stencils.
/// @param[in] ix      The grid index in x.
/// @param[in] iy      The grid index in y.
/// @param[in] iz      The grid index in z.
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridY  The number of grid points in y.
/// @param[in] nGridZ  The number of grid points in z.
/// @param[out] i0     The update node.
/// @param[out] i1     The node to the left or the right of i0.
/// @param[out] i2     The node diagonal but in the same plane as i0.
/// @param[out] i3     The node to the front of back of i0.
/// @param[out] i4     The node above or below i0.
/// @param[out] i5     The node above or below i1.
/// @param[out] i6     The node above or below i2.
/// @param[out] i7     The node above or below i3.
#pragma omp declare simd uniform(nGridX, nGridY, nGridZ)
template<EikonalXX::SweepNumber3D E>
void gridToSurroundingTravelTimeIndices(
    const int ix, const int iy, const int iz,
    const int nGridX, const int nGridY, const int nGridZ,
    int *i0, int *i1, int *i2, int *i3,
    int *i4, int *i5, int *i6, int *i7)
{
    if constexpr (E == SweepNumber3D::SWEEP1)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix - 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy - 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz - 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP2)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix + 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy - 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz - 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP3)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix - 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy + 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz - 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP4)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix + 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy + 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz - 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP5)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix - 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy - 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz + 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP6)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix + 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy - 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz + 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP7)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix - 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy + 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz + 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
    else if constexpr (E == SweepNumber3D::SWEEP8)
    {
        int ixShift = std::max(0, std::min(nGridX - 1, ix + 1));
        int iyShift = std::max(0, std::min(nGridY - 1, iy + 1));
        int izShift = std::max(0, std::min(nGridZ - 1, iz + 1));
        *i0 = gridToIndex(nGridX, nGridY, ix,      iy,      iz);
        *i1 = gridToIndex(nGridX, nGridY, ixShift, iy,      iz);
        *i2 = gridToIndex(nGridX, nGridY, ixShift, iyShift, iz);
        *i3 = gridToIndex(nGridX, nGridY, ix,      iyShift, iz);
        *i4 = gridToIndex(nGridX, nGridY, ix,      iy,      izShift);
        *i5 = gridToIndex(nGridX, nGridY, ixShift, iy,      izShift);
        *i6 = gridToIndex(nGridX, nGridY, ixShift, iyShift, izShift);
        *i7 = gridToIndex(nGridX, nGridY, ix,      iyShift, izShift);
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif 
}

/// @brief Fills the slowness vectors in each sweep with their slowness.
/// @param[in] nCellX          The number of cells in the x direction.
/// @param[in] nCellY          The number of cells in the y direction.
/// @param[in] nCellZ          The number of cells in the z direction.
/// @param[in] graph           Defines the computational graph for this sweep.
/// @param[in] slow            The cell-based slownesses in s/m.  This is an
///                            array whose dimension is [nCellX*nCellY*nCellZ].
/// @param[out] sweepSlowness  The sweep slownesses in each direction.
template<typename T, EikonalXX::SweepNumber3D E>
void slownessToSweepSlowness(const int nCellX,
                             const int nCellY,
                             const int nCellZ,
                             const Graph3D<E> &graph,
                             const T *__restrict__ slow,
                             SweepSlowness3D<T> *sweepSlowness)
{
    auto nLevels = graph.getNumberOfLevels();
    auto levelStartPtr = graph.getLevelStartPointer();
    const int *__restrict__ nodeToX = graph.getNodeInLevelToXGridPointPointer();
    const int *__restrict__ nodeToY = graph.getNodeInLevelToYGridPointPointer();
    const int *__restrict__ nodeToZ = graph.getNodeInLevelToZGridPointPointer();
    tbb::parallel_for(tbb::blocked_range<int> (0, nLevels),
                      [=](const tbb::blocked_range<int> &localLevel)
    {
        int i0, i1, ix, iy, iz;
        int iCell0, iCell1, iCell2, iCell3, iCell4, iCell5, iCell7 = 0;
        for (int level=localLevel.begin(); level != localLevel.end(); ++level)
        {
            auto i0 = levelStartPtr[level];
            auto i1 = levelStartPtr[level + 1];
            if (sweepSlowness[level].s0.size() < static_cast<size_t> (i1 - i0))
            {
                sweepSlowness[level].allocate(i1 - i0);
            }
            T *__restrict__ s0 = sweepSlowness[level].s0.data();
            T *__restrict__ s1 = sweepSlowness[level].s1.data();
            T *__restrict__ s2 = sweepSlowness[level].s2.data();
            T *__restrict__ s3 = sweepSlowness[level].s3.data();
            T *__restrict__ s4 = sweepSlowness[level].s4.data();
            T *__restrict__ s5 = sweepSlowness[level].s5.data();
            T *__restrict__ s7 = sweepSlowness[level].s7.data();
            for (int node = i0; node < i1; ++node)
            {
                ix = nodeToX[node];
                iy = nodeToY[node];
                iz = nodeToZ[node];
                gridToSurroundingSlownessIndices<E>(
                    ix, iy, iz,
                    nCellX, nCellY, nCellZ,
                    &iCell0, &iCell1, &iCell2, &iCell3,
                    &iCell4, &iCell5, &iCell7);
                s0[node - i0] = slow[iCell0];
                s1[node - i0] = slow[iCell1];
                s2[node - i0] = slow[iCell2];
                s3[node - i0] = slow[iCell3];
                s4[node - i0] = slow[iCell4];
                s5[node - i0] = slow[iCell5];
                s7[node - i0] = slow[iCell7];
            }
        }
    });
}

}
#endif
