#ifndef PRIVATE_SOLVERUTILITIES3D_HPP
#define PRIVATE_SOLVERUTILITIES3D_HPP
#include <CL/sycl.hpp>
#include <cmath>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
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
/// @brief The finite difference update for the uniform grid.
/// @param[in] factoredEikonalRadius The number of grid points in x, y, and z
///                            the update grid point must be from the source
///                            grid point to switch from the factored eikonal
///                            solver to the regular Cartesian solver.
/// @param[in] huge            The default value for an update should causality
///                            fail.
/// @param[in] sourceSlowness  The slowness at the source node in s/m.
/// @param[in] ix              The update grid point in x.
/// @param[in] iy              The update grid point in y.
/// @param[in] iz              The update grid point in z.
/// @param[in] ixShift         The grid sweeping direction shift in x when
///                            computing analytic travel times and derivatives.
/// @param[in] iyShift         The grid sweeping direction shift in y when
///                            computing analytic travel times and derivatives.
/// @param[in] izShift         The grid sweeping direction shift in x when
///                            computing analytic travel times and derivatives.
/// @param[in] iSrcX           The source's grid node in x.
/// @param[in] iSrcY           The source's grid node in y.
/// @param[in] iSrcZ           The source's grid node in z.
/// @param[in] xSourceOffset   The source's distance in x from the origin in m.
/// @param[in] ySourceOffset   The source's distance in y from the origin in m.
/// @param[in] zSourceOffset   The source's distance in z from the origin in m.
/// @param[in] s0              The slowness (s/m) in the home cell.
/// @param[in] s1              The slowness (s/m) in the cell to the left
///                            or right of the home cell.
/// @param[in] s2              The slowness (s/m) diagonal from and in the 
///                            same plane as s0.
/// @param[in] s4              The slowness (s/m) in the cell above or below
///                            the home cell. 
/// @param[in] s5              The slowness (s/m) in cell left or right of s4.
/// @param[in] s7              The slowness (s/m) in the cell in front of or
///                            in back of s4. 
/// @param[in] t1              The travel time (s) to the left or right of the
///                            update node.
/// @param[in] t2              The travel time (s) in front of or in back of the
///                            update node.
/// @param[in] t3              The travel time (s) of the node diagonal from and
///                            in the same plane as the update node.
/// @param[in] t4              The travel time (s) to the left or right of t5.
/// @param[in] t5              The travel time (s) above or below the update
///                            node.
/// @param[in] t6              The travel time (s) in front of or in back of t5.
/// @param[in] t7              The travel time (s) of the node diagonal and in
///                            the same plane as t5.
/// @result The candidate travel time at node (ix,iy,iz) in seconds.
template<typename T>
[[nodiscard]]
T finiteDifference(const int factoredEikonalRadius,
                   const T huge,
                   const T h,
                   const T sourceSlowness,
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
constexpr bool doFteik = true;
    // Default result
    T tUpd = huge;
//std::cout << ix << " " << iy << " " << iz << std::endl;
//std::cout << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << " " << t6 << " " << t7 <<std::endl;
    // Some mins I'll need for the refracted plane operators
    T minS0S3 = sycl::fmin(s0, s3); 
    T minS0S1 = sycl::fmin(s0, s1);
    T minS0S4 = sycl::fmin(s0, s4);
    // 1D operators critically refracted operators
    T t1d1 = t1 + h*sycl::fmin(minS0S3, sycl::fmin(s4, s7)); 
    T t1d2 = t3 + h*sycl::fmin(minS0S4, sycl::fmin(s1, s5));
    T t1d3 = t4 + h*sycl::fmin(minS0S1, sycl::fmin(s2, s3));
    T t1d = sycl::fmin(sycl::fmin(t1d1, t1d2), t1d3);
//std::cout << "t1d: " << t1d1 << " " << t1d2 << " " << t1d3 << std::endl;
    // Some geometric terms
    T hs03 = h*minS0S3;
    T hs01 = h*minS0S1;
    T hs04 = h*minS0S4; 
    // Cartesian
    if (std::abs(iSrcX - ix) > factoredEikonalRadius ||
        std::abs(iSrcY - iy) > factoredEikonalRadius ||
        std::abs(iSrcZ - iz) > factoredEikonalRadius)
    {
        // 2D critically refracted operators in XZ, YZ, and XY planes 
        T dtCrossXZ = t1 - t4; // 1/2*(dTdxXZ - dTzXZ)
        T detXZ = 2*(hs03*hs03) - dtCrossXZ*dtCrossXZ;
        T t2d1 = huge;
        if (detXZ > 0 && t4 <= t1 + hs03 && t1 <= t4 + hs03)
        {
            t2d1 = t5 + sycl::sqrt(detXZ);
        }
        T t2d = t2d1; //sycl::fmin(t1d, t2d1);

        T dtCrossYZ = t3 - t4; // 1/2*(dTdyYZ - dTzYZ)
        T detYZ = 2*(hs01*hs01) - dtCrossYZ*dtCrossYZ;
        T t2d2 = huge;
        if (detYZ > 0 && t3 <= t4 + hs01 && t4 <= t3 + hs01)
        {
            t2d2 = t7 + sycl::sqrt(detYZ);
        }
        t2d = sycl::fmin(t2d, t2d2);

        T dtCrossXY = t1 - t3; // 1/2*(dTdxXY - dTyXY)
        T detXY = 2*(hs04*hs04) - dtCrossXY*dtCrossXY;
        T t2d3 = huge;
        if (detXY > 0 && t3 <= t1 + hs04 && t1 <= t3 + hs04)
        {
            t2d3 = t2 + sycl::sqrt(detXY);
        }
        t2d = sycl::fmin(t2d, t2d3);
std::cout << "t2d 3point: " << t2d1 << " " << t2d2 <<  " " << t2d3 << std::endl;
//std::cout << "compare: " << t1d << " " << t2d << " " << sycl::fmax(sycl::fmax(t1, t4), t3) << std::endl;
        // 8 point operator.  Quick test to see if this is even a candidate
        // by checking if the wavefront is propagating the right way.
        T t3d = huge;
        if (sycl::fmin(t1d, t2d) > sycl::fmax(sycl::fmax(t1, t4), t3))
        {
            // This is simplified from Noble's fteik code.
            constexpr T third = 0.33333333333333333333333;
            if (doFteik)
            {
                constexpr T half = 0.5;
                T Tx = t1 + half*(-t4 + t5 - t3 + t2) - t7 + t6;
                T Ty = t3 + half*(-t4 + t7 - t1 + t2) - t5 + t6;
                T Tz = t4 + half*(-t1 + t5 - t3 + t7) - t2 + t6;
//std::cout << "ta,tb,tc: " << Tx << " " << Ty << " " << Tz << std::endl;
                T det = 27*(h*s0)*(h*s0)
                      - (Tx - Ty)*(Tx - Ty)
                      - (Tx - Tz)*(Tx - Tz)
                      - (Ty - Tz)*(Ty - Tz);
                if (det >= 0)
                {
                    t3d = third*((Tx + Ty + Tz) + sycl::sqrt(det));
                }
            }
            else
            {
                // This is the expression I could derive.  note, the
                // 48 intead of the 27 in the determinant.
                T Tx = t1 - t4 + t5 - t3 + t2 - t7 + t6;
                T Ty = t3 - t4 + t7 - t1 + t2 - t5 + t6;
                T Tz = t4 - t1 + t5 - t3 + t7 - t2 + t6;
                T det = 48*(h*s0)*(h*s0)
                      - (Tx - Ty)*(Tx - Ty)
                      - (Tx - Tz)*(Tx - Tz)
                      - (Ty - Tz)*(Ty - Tz); 
                if (det >= 0)
                {    
                    t3d = third*((Tx + Ty + Tz) + sycl::sqrt(det));
                }
            }
//std::cout << "t3d: " << t3d << std::endl;
        } 
        tUpd = sycl::fmin(t1d, sycl::fmin(t2d, t3d));
    }
    else // Spherical
    {
        T t0, dt0dx, dt0dy, dt0dz;
        computeAnalyticalTravelTime(ix, iy, iz,
                                    h, h, h,
                                    xSourceOffset, ySourceOffset, zSourceOffset,
                                    s0,
                                    &t0, &dt0dx, &dt0dy, &dt0dz);
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
 //std::cout << "tana: " << tau1 << " " << tau2 << " " << tau3 << " " << tau4 << " " << tau5 << " " << tau6 << " " << tau7 << std::endl;
        // Travel time derivative signs should be consistent with sweep
        // direction since the sweep direction is presuming the direction
        // that the wave is traveling.
        dt0dx = signX*dt0dx;
        dt0dy = signY*dt0dy; 
        dt0dz = signZ*dt0dz;
        // Analog of four point operators but for factored eikonal equation.
        T dTauTxXZ = tau1 - tau4 + tau5;
        T dTauTzXZ = tau4 - tau1 + tau5;
        //constexpr T one = 1; 
        //constexpr T two = 2; 
        //constexpr T four = 4;
        constexpr T half = 0.5;
        T hi = 1/h;
        T h2inv = 1/(h*h);
        T a = h2inv;
        T b = (2*hi)*(dt0dx + dt0dz)
            - (2*h2inv)*tau5;
        T c = (half*h2inv)*(dTauTxXZ*dTauTxXZ + dTauTzXZ*dTauTzXZ)
            - (2*hi)*(dt0dx*dTauTxXZ + dt0dz*dTauTzXZ)
            + 2*(sourceSlowness*sourceSlowness - s0*s0 + dt0dy*dt0dy);
        T detXZ = b*b - 4*a*c;
        // Critically refracted waves in XZ plane, YZ, and XY planes.
        T tau = 0;
        T t2d1 = huge;
        if (detXZ > 0 && t4 <= t1 + hs03 && t1 <= t4 + hs03)
        {
            tau = (-b + sycl::sqrt(detXZ))/(2*a);
            t2d1 = t0 + tau;
//std::cout << "potential t1: " << t2d1 << std::endl;
            // Ensure travel time is increasing 
            if (t2d1 < t1 || t2d1 < t4){t2d1 = huge;} 
        }
        T t2d = t2d1; //sycl::fmin(t1d, t2d1);

        T dTauTyYZ = tau3 - tau4 + tau7;
        T dTauTzYZ = tau4 - tau3 + tau7;
        b = (2*hi)*(dt0dy + dt0dz)
          - (2*h2inv)*tau7;
        c = (half*h2inv)*(dTauTyYZ*dTauTyYZ + dTauTzYZ*dTauTzYZ)
          - (2*hi)*(dt0dy*dTauTyYZ + dt0dz*dTauTzYZ)
          + 2*(sourceSlowness*sourceSlowness - s0*s0 + dt0dx*dt0dx);
        T detYZ = b*b - 4*a*c;
        T t2d2 = huge;
        if (detYZ >= 0 && t3 <= t4 + hs01 && t4 <= t3 + hs01)
        {
            tau = (-b + sycl::sqrt(detYZ))/(2*a);
            t2d2 = t0 + tau;
//std::cout << "potential t2: " << t2d2 << std::endl;
            if (t2d2 < t3 || t2d2 < t4){t2d2 = huge;}
        }
        t2d = sycl::fmin(t2d, t2d2);

        T dTauTxXY = tau1 - tau3 + tau2;
        T dTauTyXY = tau3 - tau1 + tau2;
        b = (2*hi)*(dt0dx + dt0dy)
          - (2*h2inv)*tau2;
        c = (half*h2inv)*(dTauTxXY*dTauTxXY + dTauTyXY*dTauTyXY)
          - (2*hi)*(dt0dx*dTauTxXY + dt0dy*dTauTyXY)
          + 2*(sourceSlowness*sourceSlowness - s0*s0 + dt0dz*dt0dz);
        T detXY = b*b - 4*a*c;
        T t2d3 = huge;
        if (detXY >= 0 && t3 <= t1 + hs04 && t1 <= t3 + hs04)
        {
            tau = (-b + sycl::sqrt(detXY))/(2*a);
            t2d3 = t0 + tau;
//std::cout << "potential t3: " << t2d3 << std::endl;
            if (t2d3 < t1 || t2d3 < t3){t2d3 = huge;}
        }
        t2d = sycl::fmin(t2d, t2d3);
//std::cout << "sph 3 point: " << t2d1 << " " << t2d2 <<  " " << t2d3 << std::endl;

        // 8 point operator
        T t3d = huge;
        if (sycl::fmin(t1d, t2d) > sycl::fmax(sycl::fmax(t1, t4), t3)) 
        {
            // This is simplified from Noble's fteik code.
            if (doFteik)
            {
                T taux = tau1 + half*(-tau4 + tau5 - tau3 + tau2) - tau7 + tau6;
                T tauy = tau3 + half*(-tau4 + tau7 - tau1 + tau2) - tau5 + tau6;
                T tauz = tau4 + half*(-tau1 + tau5 - tau3 + tau7) - tau2 + tau6;
                a = 3*h2inv;
                b =-2*h2inv*(taux + tauy + tauz)
                  + (6*hi)*(dt0dx + dt0dy + dt0dz);
                c = h2inv*(taux*taux + tauy*tauy + tauz*tauz)
                  - (6*hi)*(dt0dx*taux + dt0dy*tauy + dt0dz*tauz)
                  + 9*(sourceSlowness*sourceSlowness - s0*s0);
                T det = b*b - 4*a*c;
                if (det >= 0)
                {
                    t3d = t0 + (-b + sycl::sqrt(det))/(2*a);
                }
            }
            else
            {
                // This is the expression I derived.  Some of the coefficients
                // are slightly different.
                T taux = tau1 - tau4 + tau5 - tau3 + tau2 - tau7 + tau6;
                T tauy = tau3 - tau4 + tau7 - tau1 + tau2 - tau5 + tau6;
                T tauz = tau4 - tau1 + tau5 - tau3 + tau7 - tau2 + tau6;
                a = 3*h2inv;
                b =-2*h2inv*(taux + tauy + tauz)
                  + (8*hi)*(dt0dx + dt0dy + dt0dz);
                c = h2inv*(taux*taux + tauy*tauy + tauz*tauz)
                  - (8*hi)*(dt0dx*taux + dt0dy*tauy + dt0dz*tauz)
                  + 16*(sourceSlowness*sourceSlowness - s0*s0);
                T det = b*b - 4*a*c;
                if (det >= 0)
                {
                    t3d = t0 + (-b + sycl::sqrt(det))/(2*a);
                }
            }
            if (t3d < t1 || t3d < t3 || t3d < t4){t3d = huge;}
        }
        tUpd = sycl::fmin(t1d, sycl::fmin(t2d, t3d));
    }
std::cout << "cand tupd: (" << ix << "," << iy << "," << iz << ") " << tUpd << std::endl;
//getchar();
//if (ix == 5 && iy == 5 && iz == 3){getchar();}
    return tUpd;
}


/// @brief The finite difference for the non-uniform grid.
template<typename T>
[[nodiscard]]
T finiteDifference(const int factoredEikonalRadius,
                   const T huge,
                   const T dx, const T dy, const T dz,
                   const T dxInv, const T dyInv, const T dzInv,
                   const T dx2Inv, const T dy2Inv, const T dz2Inv,
                   const T dx_dz, const T dz_dx,
                   const T dy_dz, const T dz_dy,
                   const T dy_dx, const T dx_dy,
                   const T dx2_p_dz2_inv,
                   const T dy2_p_dz2_inv,
                   const T dx2_p_dy2_inv,
                   const T sourceSlowness,
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
constexpr bool doFteik = true;
    // Default result
    T tUpd = huge;
    // Some mins I'll need for the refracted plane operators
    T minS0S3 = sycl::fmin(s0, s3);
    T minS0S1 = sycl::fmin(s0, s1);
    T minS0S4 = sycl::fmin(s0, s4);
    // 1D operators critically refracted operators
    T t1d1 = t1 + dx*sycl::fmin(minS0S3, sycl::fmin(s4, s7));
    T t1d2 = t3 + dy*sycl::fmin(minS0S4, sycl::fmin(s1, s5));
    T t1d3 = t4 + dz*sycl::fmin(minS0S1, sycl::fmin(s2, s3));
    T t1d = sycl::fmin(sycl::fmin(t1d1, t1d2), t1d3);
    // Cartesian
    if (std::abs(iSrcX - ix) > factoredEikonalRadius ||
        std::abs(iSrcY - iy) > factoredEikonalRadius ||
        std::abs(iSrcZ - iz) > factoredEikonalRadius)
    {
        // 2D critically refracted operators in XZ, YZ, and XY planes 
        T dtCrossXZ = t1 - t4; // 1/2*(dTdxXZ - dTzXZ
        T detXZ = (dx*dx + dz*dz)*(minS0S3*minS0S3) - dtCrossXZ*dtCrossXZ;
        T t2d1 = huge;
        if (detXZ > 0 && t1 < t4 + dz*minS0S3 && t4 < t1 + dx*minS0S3)
        {
            T TxXZ = t1 - t4 + t5;
            T TzXZ = t4 - t1 + t5;
            t2d1 = (TxXZ*dz_dx + TzXZ*dx_dz) + 2*sycl::sqrt(detXZ);
            t2d1 = ((dx*dz)*dx2_p_dz2_inv)*t2d1;
        }
        T t2d = t2d1; //sycl::fmin(t1d, t2d1);

        T dtCrossYZ = t3 - t4; // 1/2*(dTdyYZ - dTzYZ)
        T detYZ = (dy*dy + dz*dz)*(minS0S1*minS0S1) - dtCrossYZ*dtCrossYZ;
        T t2d2 = huge;
        if (detYZ > 0 && t3 < t4 + dz*minS0S1 && t4 < t3 + dy*minS0S1)
        {
            T TyYZ = t3 - t4 + t7;
            T TzYZ = t4 - t3 + t7;
            t2d2 = (TyYZ*dz_dy + TzYZ*dy_dz) + 2*sycl::sqrt(detYZ);
            t2d2 = ((dy*dz)*dy2_p_dz2_inv)*t2d2;
        }
        t2d = sycl::fmin(t2d, t2d2);

        T dtCrossXY = t1 - t3; // 1/2*(dTdxXY - dTyXY)
        T detXY = (dx*dx + dy*dy)*(minS0S4*minS0S4) - dtCrossXY*dtCrossXY;
        T t2d3 = huge;
        if (detXY > 0 && t3 < t1 + dx*minS0S4 && t1 < t3 + dy*minS0S4)
        {
            T TxXY = t1 - t3 + t2;
            T TyXY = t3 - t1 + t2;
            t2d3 = (TxXY*dy_dx + TyXY*dx_dy) + 2*sycl::sqrt(detXY);
            t2d3 = ((dx*dy)*dx2_p_dy2_inv)*t2d3;
        }
        t2d = sycl::fmin(t2d, t2d3);
//std::cout << t2d1 << " " << t2d2 << " " << t2d3 << std::endl;
        // 8 point operator
        T t3d = huge;
        if (sycl::fmin(t1d, t2d) > sycl::fmax(sycl::fmax(t1, t4), t3))
        {
//std::cout << "8 point" << std::endl;
            if (doFteik)
            {
                constexpr T half = 0.5;
                T Tx = t1 + half*(-t4 + t5 - t3 + t2) - t7 + t6;
                T Ty = t3 + half*(-t4 + t7 - t1 + t2) - t5 + t6;
                T Tz = t4 + half*(-t1 + t5 - t3 + t7) - t2 + t6;
                T det = 9*(s0*s0)*(dx2Inv + dy2Inv + dz2Inv)
                      - ((Tx - Tz)*(Tx - Tz))*(dx2Inv*dz2Inv)
                      - ((Ty - Tz)*(Ty - Tz))*(dy2Inv*dz2Inv)
                      - ((Tx - Ty)*(Tx - Ty))*(dx2Inv*dy2Inv);
                if (det >= 0)
                {
                    t3d = (Tx*dx2Inv + Ty*dy2Inv + Tz*dz2Inv + sycl::sqrt(det))
                         /(dx2Inv + dy2Inv + dz2Inv);
                }
            }
            else
            {
                T Tx = t1 - t4 + t5 - t3 + t2 - t7 + t6;
                T Ty = t3 - t4 + t7 - t1 + t2 - t5 + t6;
                T Tz = t4 - t1 + t5 - t3 + t7 - t2 + t6;
                T a = dx2Inv + dy2Inv + dz2Inv;
                T b = -2*(Tx*dx2Inv + Ty*dy2Inv + Tz*dz2Inv);
                T c = (Tx*Tx)*dx2Inv + (Ty*Ty)*dy2Inv + (Tz*Tz)*dz2Inv
                    - 16*s0*s0;
                T det = b*b - 4*a*c;
                if (det >= 0)
                {
                    t3d = (-b + sycl::sqrt(det))/(2*a);
                }
            }
        }
        tUpd = sycl::fmin(t1d, sycl::fmin(t2d, t3d));
    }
    else // Spherical
    {
        T t0, dt0dx, dt0dy, dt0dz;
        computeAnalyticalTravelTime(ix, iy, iz,
                                    dx, dy, dz,
                                    xSourceOffset, ySourceOffset, zSourceOffset,
                                    s0,
                                    &t0, &dt0dx, &dt0dy, &dt0dz);
        T tau1 = t1 - computeAnalyticalTravelTime(
                          ix + ixShift, iy, iz,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau2 = t2 - computeAnalyticalTravelTime(
                          ix + ixShift, iy + iyShift, iz,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau3 = t3 - computeAnalyticalTravelTime(
                          ix, iy + iyShift, iz,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau4 = t4 - computeAnalyticalTravelTime(
                          ix, iy, iz + izShift,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau5 = t5 - computeAnalyticalTravelTime(
                          ix + ixShift, iy, iz + izShift,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau6 = t6 - computeAnalyticalTravelTime(
                          ix + ixShift, iy + iyShift, iz + izShift,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        T tau7 = t7 - computeAnalyticalTravelTime(
                          ix, iy + iyShift, iz + izShift,
                          dx, dy, dz,
                          xSourceOffset, ySourceOffset, zSourceOffset,
                          s0);
        dt0dx = signX*dt0dx;
        dt0dy = signY*dt0dy;
        dt0dz = signZ*dt0dz;

        // 4 Point Operator -> XZ
        T dTauTxXZ = tau1 - tau4 + tau5;
        T dTauTzXZ = tau4 - tau1 + tau5;
        T a = dx2Inv + dz2Inv;
        T b = 4*(dxInv*dt0dx + dzInv*dt0dz)
            - 2*(dTauTxXZ*dx2Inv + dTauTzXZ*dz2Inv);
        T c = (dTauTxXZ*dTauTxXZ)*dx2Inv + (dTauTzXZ*dTauTzXZ)*dz2Inv
            - 4*((dt0dx*dTauTxXZ)*dxInv + (dt0dz*dTauTzXZ)*dzInv)
            + 4*(sourceSlowness*sourceSlowness - s0*s0 + dt0dy*dt0dy); 
        T detXZ = b*b - 4*a*c;
        T t2d1 = huge;
        if (detXZ > 0 && t1 < t4 + dz*minS0S3 && t4 < t1 + dx*minS0S3)
        {
            t2d1 = t0 + (-b + sycl::sqrt(detXZ))/(2*a);
            if (t2d1 < t1 || t2d1 < t4){t2d1 = huge;}
        }
        T t2d = t2d1;
        // YZ
        T dTauTyYZ = tau3 - tau4 + tau7;
        T dTauTzYZ = tau4 - tau3 + tau7;
        a = dy2Inv + dz2Inv;
        b = 4*(dyInv*dt0dy + dzInv*dt0dz)
          - 2*(dTauTyYZ*dy2Inv + dTauTzYZ*dz2Inv);
        c = (dTauTyYZ*dTauTyYZ)*dy2Inv + (dTauTzYZ*dTauTzYZ)*dz2Inv
          - 4*((dt0dy*dTauTyYZ)*dyInv + (dt0dz*dTauTzYZ)*dzInv)
          + 4*(sourceSlowness*sourceSlowness - s0*s0 + dt0dx*dt0dx);
        T detYZ = b*b - 4*a*c;
        T t2d2 = huge;
        if (detYZ > 0 && t3 < t4 + dz*minS0S1 && t4 < t3 + dy*minS0S1)
        {
            t2d2 = t0 + (-b + sycl::sqrt(detYZ))/(2*a);
            if (t2d2 < t3 || t2d2 < t4){t2d2 = huge;}
        }
        t2d = sycl::fmin(t2d, t2d2);
        // XY
        T dTauTxXY = tau1 - tau3 + tau2;
        T dTauTyXY = tau3 - tau1 + tau2;
        a = dx2Inv + dy2Inv;
        b = 4*(dxInv*dt0dx + dyInv*dt0dy)
          - 2*(dTauTxXY*dx2Inv + dTauTyXY*dy2Inv);
        c = (dTauTxXY*dTauTxXY)*dx2Inv + (dTauTyXY*dTauTyXY)*dy2Inv
          - 4*((dt0dx*dTauTxXY)*dxInv + (dt0dy*dTauTyXY)*dyInv)
          + 4*(sourceSlowness*sourceSlowness - s0*s0 + dt0dz*dt0dz); 
        T detXY = b*b - 4*a*c;
        T t2d3 = huge;
        if (detXY > 0 && t3 < t1 + dx*minS0S4 && t1 < t3 + dy*minS0S4)
        {
            t2d3 = t0 + (-b + sycl::sqrt(detXY))/(2*a);
            if (t2d3 < t1 || t2d3 < t3){t2d3 = huge;}
        }
        t2d = sycl::fmin(t2d, t2d3);
        // 8 point operator
        T t3d = huge;
        if (sycl::fmin(t1d, t2d) > sycl::fmax(sycl::fmax(t1, t4), t3))
        {
            if (doFteik)
            {
                constexpr T half = 0.5;
                T taux = tau1 + half*(-tau4 + tau5 - tau3 + tau2) - tau7 + tau6;
                T tauy = tau3 + half*(-tau4 + tau7 - tau1 + tau2) - tau5 + tau6;
                T tauz = tau4 + half*(-tau1 + tau5 - tau3 + tau7) - tau2 + tau6;
                a = dx2Inv + dy2Inv + dz2Inv;
                b =-2*(taux*dx2Inv + tauy*dy2Inv + tauz*dz2Inv)
                  + 6*(dxInv*dt0dx + dyInv*dt0dy + dzInv*dt0dz);
                c = (taux*taux)*dx2Inv + (tauy*tauy)*dy2Inv + (tauz*tauz)*dz2Inv
                  - 6*(dt0dx*taux*dxInv + dt0dy*tauy*dyInv + dt0dz*tauz*dzInv)
                  + 9*(sourceSlowness*sourceSlowness - s0*s0);
                T det = b*b - 4*a*c;
                if (det >= 0)
                {
                    t3d = t0 + (-b + sycl::sqrt(det))/(2*a);
                }
            }
            else
            {
                T taux = tau1 - tau4 + tau5 - tau3 + tau2 - tau7 + tau6;
                T tauy = tau3 - tau4 + tau7 - tau1 + tau2 - tau5 + tau6;
                T tauz = tau4 - tau1 + tau5 - tau3 + tau7 - tau2 + tau6;
                a = dx2Inv + dy2Inv + dz2Inv;
                b =-2*(taux*dx2Inv + tauy*dy2Inv + tauz*dz2Inv)
                  + 8*(dxInv*dt0dx + dyInv*dt0dy + dzInv*dt0dz);
                c = (taux*taux)*dx2Inv + (tauy*tauy)*dy2Inv + (tauz*tauz)*dz2Inv
                  - 8*(dt0dx*taux*dxInv + dt0dy*tauy*dyInv + dt0dz*tauz*dzInv)
                  + 16*(sourceSlowness*sourceSlowness - s0*s0);
                T det = b*b - 4*a*c;
                if (det >= 0)
                {
                    t3d = t0 + (-b + sycl::sqrt(det))/(2*a);
                }
            }
            if (t3d < t1 || t3d < t3 || t3d < t4){t3d = huge;}
        }
        tUpd = sycl::fmin(t1d, sycl::fmin(t2d, t3d));
    }
std::cout << "anisotropic cand tupd: (" << ix << "," << iy << "," << iz << ") " << tUpd << std::endl;
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
    if constexpr(E == EikonalXX::SweepNumber3D::Sweep1)
    {
        *ixShift =-1;
        *iyShift =-1;
        *izShift =-1;
        *signX = 1;
        *signY = 1;
        *signZ = 1; 
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep2)
    {
        *ixShift =+1;
        *iyShift =-1;
        *izShift =-1;
        *signX =-1;
        *signY = 1;
        *signZ = 1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep3)
    {
        *ixShift =-1;
        *iyShift =+1;
        *izShift =-1;
        *signX = 1;
        *signY =-1;
        *signZ = 1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep4)
    {
        *ixShift =+1;
        *iyShift =+1;
        *izShift =-1;
        *signX =-1;
        *signY =-1;
        *signZ = 1; 
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep5) // Sweep 1; flip z
    {
        *ixShift =-1;
        *iyShift =-1;
        *izShift =+1;
        *signX = 1;
        *signY = 1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep6) // Sweep 2; flip z
    {
        *ixShift =+1;
        *iyShift =-1;
        *izShift =+1;
        *signX =-1;
        *signY = 1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep7) // Sweep 3; flip z
    {
        *ixShift =-1;
        *iyShift =+1;
        *izShift =+1;
        *signX = 1;
        *signY =-1;
        *signZ =-1;
    }
    else if constexpr(E == EikonalXX::SweepNumber3D::Sweep8) // Sweep 4 ; flip z
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
    if constexpr (E == SweepNumber3D::Sweep1)
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
    else if constexpr (E == SweepNumber3D::Sweep2)
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
    else if constexpr (E == SweepNumber3D::Sweep3)
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
    else if constexpr (E == SweepNumber3D::Sweep4)
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
    else if constexpr (E == SweepNumber3D::Sweep5)
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
    else if constexpr (E == SweepNumber3D::Sweep6)
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
    else if constexpr (E == SweepNumber3D::Sweep7)
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
    else if constexpr (E == SweepNumber3D::Sweep8)
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

/// @brief Gets the loop limits for the initialization sweep of the
///        fast sweeping method i.e.,:
///        for (int iz = iz0; iz != iz1; iz = iz + izDir)
///            for (int iy = iy0; iy != iy1; iy = iy + iyDir)
///                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
/// @param[in] sweep   The sweep number.
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridY  The number of grid points in y.
/// @param[in] nGridZ  The number of grid points in z.
/// @param[in] iSrcX   The source grid point in x.
/// @param[in] iSrcY   The source grid point in y.
/// @param[in] iSrcZ   The source grid point in z.
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
                   const int iSrcX, const int iSrcY, const int iSrcZ, 
                   int *ix0, int *iy0, int *iz0,
                   int *ix1, int *iy1, int *iz1,
                   int *ixDir, int *iyDir, int *izDir)
{
    if constexpr (E == SweepNumber3D::Sweep1)
    {
        *ix0 = std::max(1, iSrcX);
        *ix1 = nGridX;
        *iy0 = std::max(1, iSrcY);
        *iy1 = nGridY;
        *iz0 = std::max(1, iSrcZ);
        *iz1 = nGridZ;
        *ixDir = 1;
        *iyDir = 1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::Sweep2)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iy0 = std::max(1, iSrcY);
        *iy1 = nGridY;
        *iz0 = std::max(1, iSrcZ);
        *iz1 = nGridZ;
        *ixDir =-1;
        *iyDir = 1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::Sweep3)
    {
        *ix0 = std::max(1, iSrcX);
        *ix1 = nGridX;
        *iy0 = std::min(iSrcY + 1, nGridY - 2);
        *iy1 =-1;
        *iz0 = std::max(1, iSrcZ);
        *iz1 = nGridZ;
        *ixDir = 1;
        *iyDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::Sweep4)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iy0 = std::min(iSrcY + 1, nGridY - 2);
        *iy1 =-1;
        *iz0 = std::max(1, iSrcZ);
        *iz1 = nGridZ;
        *ixDir =-1;
        *iyDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber3D::Sweep5)
    {
        *ix0 = std::max(1, iSrcX);
        *ix1 = nGridX;
        *iy0 = std::max(1, iSrcY);
        *iy1 = nGridY;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
        *iz1 =-1;
        *ixDir = 1;
        *iyDir = 1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::Sweep6)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iy0 = std::max(1, iSrcY);
        *iy1 = nGridY;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
        *iz1 =-1;
        *ixDir =-1;
        *iyDir = 1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::Sweep7)
    {
        *ix0 = std::max(1, iSrcX);
        *ix1 = nGridX;
        *iy0 = std::min(iSrcY + 1, nGridY - 2);
        *iy1 =-1;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
        *iz1 =-1;
        *ixDir = 1;
        *iyDir =-1;
        *izDir =-1;
    }
    else if constexpr (E == SweepNumber3D::Sweep8)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iy0 = std::min(iSrcY + 1, nGridY - 2);
        *iy1 =-1;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
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
    if constexpr (E == SweepNumber3D::Sweep1)
    {
        *jx = ix;
        *jy = iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep2)
    {
        *jx = nx - 1 - ix;
        *jy = iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep3)
    {
        *jx = ix;
        *jy = ny - 1 - iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep4)
    {
        *jx = nx - 1 - ix;
        *jy = ny - 1 - iy;
        *jz = iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep5)
    {
        *jx = ix;
        *jy = iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep6)
    {
        *jx = nx - 1 - ix;
        *jy = iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep7)
    {
        *jx = ix;
        *jy = ny - 1 - iy;
        *jz = nz - 1 - iz;
    }
    else if constexpr (E == SweepNumber3D::Sweep8)
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
    if constexpr (E == SweepNumber3D::Sweep1)
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
    else if constexpr (E == SweepNumber3D::Sweep2)
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
    else if constexpr (E == SweepNumber3D::Sweep3)
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
    else if constexpr (E == SweepNumber3D::Sweep4)
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
    else if constexpr (E == SweepNumber3D::Sweep5)
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
    else if constexpr (E == SweepNumber3D::Sweep6)
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
    else if constexpr (E == SweepNumber3D::Sweep7)
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
    else if constexpr (E == SweepNumber3D::Sweep8)
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
    if constexpr (E == SweepNumber3D::Sweep1)
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
    else if constexpr (E == SweepNumber3D::Sweep2)
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
    else if constexpr (E == SweepNumber3D::Sweep3)
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
    else if constexpr (E == SweepNumber3D::Sweep4)
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
    else if constexpr (E == SweepNumber3D::Sweep5)
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
    else if constexpr (E == SweepNumber3D::Sweep6)
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
    else if constexpr (E == SweepNumber3D::Sweep7)
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
    else if constexpr (E == SweepNumber3D::Sweep8)
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

/// @brief Determines if the grid point is on the sweep's boundary - i.e., 
///        this grid point should not be updated.
/// @param[in] ix   The grid point in x.
/// @param[in] iy   The grid point in y.
/// @param[in] iz   The grid point in z.
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridY  The number of grid points in y.
/// @param[in] nGridZ  The number of grid points in z.
/// @result True indicates that the (ix,iy,iz) grid points is on the
///         sweep's boundary.
#pragma omp declare simd uniform(nGridX, nGridY, nGridZ)
template<EikonalXX::SweepNumber3D E>
bool isSweepBoundaryNode(const int ix, const int iy, const int iz,
                         const int nGridX, const int nGridY, const int nGridZ)
{
    bool onBoundary = false;
    if constexpr (E == SweepNumber3D::Sweep1)
    {
        onBoundary = (ix == 0 || iy == 0 || iz == 0);
    }
    else if constexpr (E == SweepNumber3D::Sweep2)
    {
        onBoundary = (ix == nGridX - 1 || iy == 0 || iz == 0);
    }
    else if constexpr (E == SweepNumber3D::Sweep3)
    {
        onBoundary = (ix == 0 || iy == nGridY - 1 || iz == 0);
    }
    else if constexpr (E == SweepNumber3D::Sweep4)
    {
        onBoundary = (ix == nGridX - 1 || iy == nGridY - 1 || iz == 0);
    }
    else if constexpr (E == SweepNumber3D::Sweep5)
    {
        onBoundary = (ix == 0 || iy == 0 || iz == nGridZ - 1);
    }
    else if constexpr (E == SweepNumber3D::Sweep6)
    {
        onBoundary = (ix == nGridX - 1 || iy == 0 || iz == nGridZ - 1);
    }
    else if constexpr (E == SweepNumber3D::Sweep7)
    {
        onBoundary = (ix == 0 || iy == nGridY - 1 || iz == nGridZ - 1);
    }
    else if constexpr (E == SweepNumber3D::Sweep8)
    {
        onBoundary = (ix == nGridX - 1 || iy == nGridY - 1 || iz == nGridZ - 1);
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
    return onBoundary;
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
            i0 = levelStartPtr[level];
            i1 = levelStartPtr[level + 1];
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
