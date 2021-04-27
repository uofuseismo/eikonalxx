#ifndef PRIVATE_SOLVERUTILITIES2D_HPP
#define PRIVATE_SOLVERUTILITIES2D_HPP
#include <CL/sycl.hpp>
#include <cmath>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#include <tbb/tbb.h>
#include "private/grid.hpp"
#include "private/pad.hpp"

#define BOUNDARY_NODE -1
#define SOURCE_NODE 0
#define UPDATE_NODE 1
namespace
{

using namespace EikonalXX;
/// Defines the sweep number
/*
enum class SweepNumber : int
{
    Sweep1 = 0,
    Sweep2 = 1,
    Sweep3 = 2,
    Sweep4 = 3
};
*/

/// @brief Container for the slownesses in a level for the 2D solver.
template<class T>
struct SweepSlowness2D
{
    /// C'tor
    SweepSlowness2D() = default;
    /// Copy c'tor
    SweepSlowness2D(const SweepSlowness2D &slowness)
    {
        *this = slowness;
    }
    /// Move c'tor
    SweepSlowness2D(SweepSlowness2D &&slowness) noexcept
    {
        *this = std::move(slowness);
    }
    /// Copy assingment
    SweepSlowness2D& operator=(const SweepSlowness2D &slowness)
    {
        if (&slowness == this){return *this;}
        if (slowness.nNodes > 0)
        {
            allocate(nNodes);
            std::copy(slowness.s0, s0, nNodes);
            std::copy(slowness.s1, s1, nNodes);
            std::copy(slowness.s3, s3, nNodes);
        }
        else
        {
            nNodes = 0;
            s0 = nullptr;
            s1 = nullptr;
            s3 = nullptr;
        }
        return *this;
    }
    /// Move assignment
    SweepSlowness2D& operator=(SweepSlowness2D &&slowness) noexcept
    {
        if (&slowness == this){return *this;}
        nNodes = slowness.nNodes; 
        s0 = std::move(slowness.s0); 
        s1 = std::move(slowness.s1);
        s3 = std::move(slowness.s3);
        return *this;
    }
    /// Destructor
    ~SweepSlowness2D()
    {
        clear();
    }
    /// Allocates space to hold slownesses for the n nodes in the sweep.
    void allocate(int n)
    {
        clear();
        if (n <= 0){throw std::invalid_argument("n must be positive");}
        nNodes = n;
        auto nPad = padLength(n, sizeof(T), 64);
        auto nBytes = nPad*sizeof(T);
        s0 = static_cast<T *> (std::aligned_alloc(64, nBytes)); 
        s1 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s3 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        std::fill(s0, s0 + nPad, 0);
        std::fill(s1, s1 + nPad, 0);
        std::fill(s3, s3 + nPad, 0);
    }
    /// This sets all the slownesses in the level to 0
    void zero() noexcept
    {
        if (nNodes > 0)
        {
            std::fill(s0, s0 + nNodes, 0);
            std::fill(s1, s1 + nNodes, 0);
            std::fill(s3, s3 + nNodes, 0);
        }
    }
    /// Gets the min slowness in the sweep
    T getMinimumValue() const noexcept    
    {
        T sMin = 0;
        if (nNodes > 0)
        {
            sMin = *std::min_element(s0, s0 + nNodes);
            sMin = std::min(sMin, *std::min_element(s1, s1 + nNodes));
            sMin = std::min(sMin, *std::min_element(s3, s3 + nNodes));
        }
        return sMin;
    } 
    /// Releases memory.
    void clear() noexcept
    {
        if (s0){free(s0);}
        if (s1){free(s1);}
        if (s3){free(s3);}
        s0 = nullptr;
        s1 = nullptr;
        s3 = nullptr;
        nNodes = 0;
    }
    /// The slowness (s/m) in the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s0 = nullptr;
    /// The slowness (s/m) in the cell to the right or left of the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s1 = nullptr;
    /// The slowness (s/m) in the cell above or below the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s3 = nullptr;
    /// The number of nodes in the sweep.
    int nNodes = 0;
};

///--------------------------------------------------------------------------///
///                             Analytical Solutions                         ///
///--------------------------------------------------------------------------///
/// @brief Computes the travel time from a source to a grid point in a 
///        homogeneous media.
/// @param[in] ix             The x grid index.
/// @param[in] iz             The z grid index.
/// @param[in] dx             The grid spacing in x in meters.
/// @param[in] dz             The grid spacing in z in meters.
/// @param[in] xSourceOffset  The source's x offset from the origin in meters.
///                           This is the source x location - model x origin.
/// @param[in] zSourceOffset  The source's z offset from the origin in meters.
///                           This is the source z location - model z origin.
/// @param[in] slowness       The slowness at the source in s/m. 
/// @result The travel time from the source to the grid point in seconds.
#pragma omp declare simd uniform(dx, dz, xSourceOffset, zSourceOffset) 
template<class T>
[[nodiscard]]
T computeAnalyticalTravelTime(const int ix, const int iz,
                              const T dx, const T dz,
                              const T xSourceOffset,
                              const T zSourceOffset,
                              const T slowness)
{
    auto x = ix*dx;
    auto z = iz*dz;
    auto deltaX = x - xSourceOffset;
    auto deltaZ = z - zSourceOffset;
    auto distance = sycl::hypot(deltaX, deltaZ); // sqrt(x**2 + z**2)
    return slowness*distance;
}
/// @brief Computes the travel time from a source to a grid point at as
///        well as the derivatives at the grid point in a homogeneous media.
/// @param[in] ix             The x grid index.
/// @param[in] iz             The z grid index.
/// @param[in] dx             The grid spacing in x in meters.
/// @param[in] dz             The grid spacing in z in meters.
/// @param[in] xSourceOffset  The source's x offset from the origin in meters.
///                           This is the source x location - model x origin.
/// @param[in] zSourceOffset  The source's z offset from the origin in meters.
///                           This is the source z location - model z origin.
/// @param[in] slowness       The slowness at the source in s/m.
/// @param[out] t             The travel time from the source to the grid
///                           point in seconds.
/// @param[out] dtdx          The derivative of the travel time w.r.t. x
///                           at the grid point in s/m. 
/// @param[out] dtdz          The derivative of the travel time w.r.t. z
///                           at the grid point in s/m.
#pragma omp declare simd uniform(dx, dz, xSourceOffset, zSourceOffset)
template<class T>
void computeAnalyticalTravelTime(const int ix, const int iz,
                                 const T dx, const T dz,
                                 const T xSourceOffset,
                                 const T zSourceOffset,
                                 const T slowness,
                                 T *t, T *dtdx, T *dtdz)
{
    auto x = ix*dx;
    auto z = iz*dz;
    auto deltaX = x - xSourceOffset;
    auto deltaZ = z - zSourceOffset;
    auto distance = sycl::hypot(deltaX, deltaZ); // sqrt(x**2 + z**2)
    *t = slowness*distance;
    // Limit x -> 0 from the right of x/sqrt(x) is 0.
    *dtdx = 0; //std::numeric_limits<T>::quiet_NaN();
    *dtdz = 0; //std::numeric_limits<T>::quiet_NaN();
    if (distance > 0)
    {
        *dtdx = (slowness*deltaX)/distance;
        *dtdz = (slowness*deltaZ)/distance; 
    }
}

///--------------------------------------------------------------------------///
///                        Cartesian Finite Differences                      ///
///--------------------------------------------------------------------------///

/*
template<typename T>
[[nodiscard]] 
T cartesianFiniteDifference(const T huge,
                            const T dx, const T dz,
                            const T dx_dz, const T dz_dx,
                            const T cosTheta, //dx/sqrt(dx*dx + dz*dz)
                            const T sinTheta, //dz/sqrt(dx*dx + dz*dz)
                            const T s0, const T s1, const T s3,
                            const T t1, const T t2, const T t3)
{
    // Critically refracted wave operators
    T t0x = t1 + dx*sycl::fmin(s0, s3);
    T t0z = t3 + dz*sycl::fmin(s0, s1);
    T t0Update1 = sycl::fmin(t0x, t0z);

    // Two three-point operators (Podvin and Lecomte, 1991)
    T hzs0 = dz*s0;
    T dTz = t1 - t2;
    T sz = s0*sinTheta;
    T t0a = huge;
    T deta = hzs0*hzs0 - dTz*dTz;
    if (dTz >= 0 && deta >= 0 && dTz <= dz*sz)
    {
        t0a = t1 + dx_dz*sycl::sqrt(deta);
    }

    T hxs0 = dx*s0;
    T dTx = t3 - t2;
    T sx = s0*cosTheta;
    T t0b = huge;
    T detb = hxs0*hxs0 - dTx*dTx;
    if (dTx >= 0 && detb >= 0 && dTx <= dx*sx)
    {
        t0b = t3 + dz_dx*sycl::sqrt(detb);
    }
    T t0Update2 = sycl::fmin(t0a, t0b);

    // Three-point operator (Vidale, 1988) 
    // Note the following identity
    // (dx*dz)/(dx^2 + dz^2) = dx/sqrt(dx^2 + dz^2) * dz/sqrt(dx^2 + dz^2)
    T dxdz_dx2_p_dz2 = cosTheta*sinTheta;
    T tildeTx = t1 - t3 + t2;
    T tildeTz = t3 - t1 + t2;
    T dtCross = t1 - t3; // tildeX - tildeZ = 2*(t1 - t3)
    T detc = hxs0*hxs0 + hzs0*hzs0 - dtCross*dtCross;
    T t0Update3 = huge;
    if (dTx >= 0 && dTz >= 0 && detc >= 0 && t1 < t3 + dz*s0 && t3 < t1 + dx*s0)
    {
        // Note for the first term:
        // (tx*dz^2 + tz*dx^2)/(dx^2 + dz^2)
        //=((tx*dz/dx + tz*dx/dz)*(dz*dx))/(dx^2 + dz^2)
        //=(tx*dz/dx + tz*dx/dz)*[(dx*dz)/(dx^2 + dz^2)]
        // Therefore, the (dx*dz)/(dx^2 + dz^2) can be factored out.
        t0Update3 = tildeTx*dz_dx + tildeTz*dx_dz + 2*sycl::sqrt(detc);
        t0Update3 = dxdz_dx2_p_dz2*t0Update3;
//std::cout << "t0Update3: " << t1 << " " << t2 << " " << t3 << " " << t0Update3  << std::endl;//" " << dxdz_dx2_p_dz2 << std::endl;//(tildeTx*dz_dx + tildeTz*dx_dz)*dxdz_dx2_p_dz2 << std::endl;
    }
//if (t0Update2 < t0Update3){std::cout<< " " << t0Update2 << " " << t0Update3 << " " <<std::endl; getchar();}
    if (t0Update3 < 0){t0Update3 = huge;}
    // Return smallest of update travel times
    return sycl::fmin(sycl::fmin(t0Update1, t0Update2), t0Update3);
}

template<typename T>
[[nodiscard]]
T cartesianFiniteDifference(const T huge,
                            const T h,
                            const T s0, const T s1, const T s3,
                            const T t1, const T t2, const T t3)
{
    // Critically refracted wave operators
    T t0x = t1 + h*sycl::fmin(s0, s3);
    T t0z = t3 + h*sycl::fmin(s0, s1);
    T t0Update1 = sycl::fmin(t0x, t0z);

    // Two-point operators (Podvin and Lecomte, 1991)
    T hs0 = h*s0;
    T dTz = t1 - t2;
    T sz = M_SQRT1_2*s0;
    T t0a = huge;
    T deta = hs0*hs0 - dTz*dTz;
    if (dTz >= 0 && deta >= 0 && dTz <= h*sz)
    {
        t0a = t1 + h*sycl::sqrt(deta);
    }

    T dTx = t3 - t2;
    T sx = M_SQRT1_2*s0;
    T t0b = huge;
    T detb = hs0*hs0 - dTx*dTx;
    if (dTx >= 0 && detb >= 0 && dTx <= h*sx)
    {
        t0b = t3 + h*sycl::sqrt(detb);
    }
    T t0Update2 = sycl::fmin(t0a, t0b);

    // Three-point operator (Vidale, 1988)
    T tildeTx = t1 - t3 + t2;
    T tildeTz = t3 - t1 + t2;
    T dtCross = t1 - t3; // tildeX - tildeZ = 2*t1 - t3
    T detc = 2*(hs0*hs0) - dtCross*dtCross;
    T t0Update3 = huge;
    if (dTx >= 0 && dTz >= 0 && detc >= 0 && t3 <= t1 + hs0 && t1 <= t3 + hs0)
    {
        t0Update3 = t2 + sycl::sqrt(detc); 
    }
    //else
    //{
    //    t0Update2 = huge;
    //    if (dTz > 0 && deta >= 0 && dTz <= h*sz)
    //    {

    //        t0Update3 = t1 + h*sycl::sqrt(deta); 
    //    }
    //    else if (dTx > 0 && detb >= 0 && dTx <= h*sx)
    //    {
    //        t0Update3 = t3 + h*sycl::sqrt(detb);
    //    }
    //}
    return sycl::fmin(sycl::fmin(t0Update1, t0Update2), t0Update3);
}
*/

template<typename T>
[[nodiscard]]
T finiteDifference(const int factoredEikonalRadius,
                   const T huge,
                   const T h,
                   const T sourceSlowness,
                   const int ix, const int iz,
                   const int signX, const int signZ,
                   const int ixShift, const int izShift, 
                   const int iSrcX, const int iSrcZ,
                   const T xSourceOffset, const T zSourceOffset,
                   const T s0, const T s1, const T s3,
                   const T t1, const T t2, const T t3)
{
    // Default result
    T tUpd = huge;
    // Critically refracted wave operators
    T t0x = t1 + h*sycl::fmin(s0, s3);
    T t0z = t3 + h*sycl::fmin(s0, s1);
    T t0Update1 = sycl::fmin(t0x, t0z);
    // Things I'll need either way
    T hs0 = h*s0;
    T dTx = t3 - t2;
    T dTz = t1 - t2;
    T sx = M_SQRT1_2*s0;
    T sz = M_SQRT1_2*s0;
    // Cartesian
    if (std::abs(iSrcX - ix) > factoredEikonalRadius ||
        std::abs(iSrcZ - iz) > factoredEikonalRadius)
    {
        // Two-point operators (Podvin and Lecomte, 1991)
        T t0a = huge;
        T deta = hs0*hs0 - dTz*dTz;
        if (dTz >= 0 && deta >= 0 && dTz <= h*sz)
        {
            t0a = t1 + h*sycl::sqrt(deta);
        }

        T t0b = huge;
        T detb = hs0*hs0 - dTx*dTx;
        if (dTx >= 0 && detb >= 0 && dTx <= h*sx)
        {
            t0b = t3 + h*sycl::sqrt(detb);
        }
        T t0Update2 = sycl::fmin(t0a, t0b);

        // Three-point operator (Vidale, 1988).  Note, we require detc > 0
        // since detc = 0 -> t0Update3 = t2 which will not result in an update.
        T dtCross = t1 - t3; // 1/2*(tildeX - tildeZ)
        T detc = 2*(hs0*hs0) - dtCross*dtCross;
        T t0Update3 = huge;
        if (dTx >= 0 && dTz >= 0 && detc > 0 &&
            t3 <= t1 + hs0 && t1 <= t3 + hs0)
        {
            t0Update3 = t2 + sycl::sqrt(detc);
        }
        tUpd = sycl::fmin(sycl::fmin(t0Update1, t0Update2), t0Update3);
    }
    else // Spherical update
    {
        T t0Update3 = huge;
        // Verify wave is moving in right direction and can reach
        //if (t1 >= t2 && t3 >= t2 && t3 < t1 + h*sx && t1 < t3 + h*sz)
        if (dTx >= 0 && dTz >= 0 && t1 < t3 + hs0 && t3 < t1 + hs0)
        {
            T dt0dx, dt0dz, t0;
            computeAnalyticalTravelTime(ix, iz, h, h,
                                        xSourceOffset, zSourceOffset,
                                        sourceSlowness,
                                        &t0, &dt0dx, &dt0dz);
            // Compute perturbations t = t_0 + tau -> tau = t - t_0
            T tau1 = t1 - computeAnalyticalTravelTime(ix + ixShift,
                                                      iz,
                                                      h, h,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness); 
            T tau2 = t2 - computeAnalyticalTravelTime(ix + ixShift,
                                                      iz + izShift,
                                                      h, h,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness);
            T tau3 = t3 - computeAnalyticalTravelTime(ix,
                                                      iz + izShift,
                                                      h, h,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness);
//std::cout << "tau: " << tau1 <<  " " << tau2 << " " << tau3  << " " << xSourceOffset << " " << zSourceOffset << std::endl;
//std::cout << "t: " << t1 << " " << t2 << " " << t3 << std::endl;
//std::cout << "twork: " << twork1 << " " << twork2 << " " << twork3 << std::endl;
            // Travel time derivative signs should be consistent with sweep
            // direction since the sweep direction is presuming the direction
            // that the wave is traveling.
            dt0dx = signX*dt0dx;
            dt0dz = signZ*dt0dz;
            // Analog of four point operators
            T tauX = tau1 - tau3 + tau2;
            T tauZ = tau3 - tau1 + tau2;
//std::cout << "taux,tauz: " << dt0dx << " " << dt0dz << std::endl;
            constexpr T one = 1;
            constexpr T two = 2;
            constexpr T four = 4;
            T h2 = h*h;
            T a = one/h2;
            T b = two/h*(dt0dx + dt0dz) - two/h2*tau2;
            T c = one/(two*h2)*(tauX*tauX + tauZ*tauZ)
                - two/h*(dt0dx*tauX + dt0dz*tauZ)
                + two*(sourceSlowness*sourceSlowness - s0*s0);
            T det = b*b - four*a*c;
            if (det >= 0)
            {
                T tau = (-b + sycl::sqrt(det))/(two*a);
//std::cout << t0 << std::endl;
//std::cout << tau << std::endl;
//std::cout << tau + t0 << std::endl;
                t0Update3 = t0 + tau;
            }
            // Ensure travel time is increasing
            if (t0Update3 < t1 || t0Update3 < t3){t0Update3 = huge;}
//std::cout<< "t0Update3 " << t0Update3 << std::endl;
        }
        tUpd = sycl::fmin(t0Update1, t0Update3);
    }
//std::cout << "t_new=" << tUpd << std::endl;
    return tUpd;
}

template<typename T>
[[nodiscard]]
T finiteDifference(const int factoredEikonalRadius,
                   const T huge,
                   const T dx, const T dz,
                   const T dx_dz, const T dz_dx,
                   const T cosTheta, const T sinTheta,
                   const T sourceSlowness,
                   const int ix, const int iz,
                   const int signX, const int signZ,
                   const int ixShift, const int izShift, 
                   const int iSrcX, const int iSrcZ,
                   const T xSourceOffset, const T zSourceOffset,
                   const T s0, const T s1, const T s3,
                   const T t1, const T t2, const T t3)
{
    T tUpd = huge;
    // Critically refracted wave operators
    T t0x = t1 + dx*sycl::fmin(s0, s3); 
    T t0z = t3 + dz*sycl::fmin(s0, s1); 
    T t0Update1 = sycl::fmin(t0x, t0z);
    // Things I'll need either way
    T dTx = t3 - t2;
    T dTz = t1 - t2;
    T sx = cosTheta*s0;
    T sz = sinTheta*s0;
    T hxs0 = dx*s0;
    T hzs0 = dz*s0;
    // Cartesian
    if (std::abs(iSrcX - ix) > factoredEikonalRadius ||
        std::abs(iSrcZ - iz) > factoredEikonalRadius)
    {
        // Two three-point operators (Podvin and Lecomte, 1991)
        T t0a = huge;
        T deta = hzs0*hzs0 - dTz*dTz;
        if (dTz >= 0 && deta >= 0 && dTz <= dz*sz)
        {
            t0a = t1 + dx_dz*sycl::sqrt(deta);
        }

        T t0b = huge;
        T detb = hxs0*hxs0 - dTx*dTx;
        if (dTx >= 0 && detb >= 0 && dTx <= dx*sx)
        {
            t0b = t3 + dz_dx*sycl::sqrt(detb);
        } 
        T t0Update2 = sycl::fmin(t0a, t0b);

        // Three-point operator (Vidale, 1988).  As before, we require detc > 0.
        // Note the following identity
        // (dx*dz)/(dx^2 + dz^2) = dx/sqrt(dx^2 + dz^2) * dz/sqrt(dx^2 + dz^2)
        T dxdz_dx2_p_dz2 = cosTheta*sinTheta;
        T tildeTx = t1 - t3 + t2;
        T tildeTz = t3 - t1 + t2;
        T dtCross = t1 - t3; // 1/2*(tildeX - tildeZ)
        T detc = hxs0*hxs0 + hzs0*hzs0 - dtCross*dtCross;
        T t0Update3 = huge; 
        if (dTx >= 0 && dTz >= 0 && detc > 0 &&
            t1 < t3 + dz*s0 && t3 < t1 + dx*s0)
        {
            // Note for the first term:
            // (tx*dz^2 + tz*dx^2)/(dx^2 + dz^2)
            //=((tx*dz/dx + tz*dx/dz)*(dz*dx))/(dx^2 + dz^2)
            //=(tx*dz/dx + tz*dx/dz)*[(dx*dz)/(dx^2 + dz^2)]
            // Therefore, the (dx*dz)/(dx^2 + dz^2) can be factored out.
            t0Update3 = tildeTx*dz_dx + tildeTz*dx_dz + 2*sycl::sqrt(detc);
            t0Update3 = dxdz_dx2_p_dz2*t0Update3;
//std::cout << "t0Update3: " << t1 << " " << t2 << " " << t3 << " " << t0Update3  << std::endl;//" " << dxdz_dx2_p_dz2 << std::endl;//(tildeTx*dz_dx + tildeTz*dx_dz)*dxdz_dx2_p_dz2 << std::endl;
        }
//if (t0Update2 < t0Update3){std::cout<< " " << t0Update2 << " " << t0Update3 << " " <<std::endl; getchar();}
        if (t0Update3 < 0){t0Update3 = huge;}
        // Return smallest of update travel times
        tUpd = sycl::fmin(sycl::fmin(t0Update1, t0Update2), t0Update3);
    }
    else // Spherical
    {
        T t0Update3 = huge;
        // Verify wave is moving in right direction and can reach
        if (dTx >= 0 && dTz >= 0 && t1 < t3 + hxs0 && t3 < t1 + hzs0) 
        {
            T dt0dx, dt0dz, t0;
            computeAnalyticalTravelTime(ix, iz, dx, dz,
                                        xSourceOffset, zSourceOffset,
                                        sourceSlowness,
                                        &t0, &dt0dx, &dt0dz);
            // Compute perturbations t = t_0 + tau -> tau = t - t_0
            T tau1 = t1 - computeAnalyticalTravelTime(ix + ixShift,
                                                      iz,
                                                      dx, dz,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness);
            T tau2 = t2 - computeAnalyticalTravelTime(ix + ixShift,
                                                      iz + izShift,
                                                      dx, dz,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness);
            T tau3 = t3 - computeAnalyticalTravelTime(ix,
                                                      iz + izShift,
                                                      dx, dz,
                                                      xSourceOffset,
                                                      zSourceOffset,
                                                      sourceSlowness);
            // Travel time derivative signs should be consistent with sweep
            // direction since the sweep direction is presuming the direction
            // that the wave is traveling.
            dt0dx = signX*dt0dx;
            dt0dz = signZ*dt0dz;
            // Analog of four point operators
            T tauX = tau1 - tau3 + tau2;
            T tauZ = tau3 - tau1 + tau2;

            constexpr T one = 1;
            constexpr T two = 2; 
            constexpr T four = 4; 
            T hx2inv = one/(dx*dx);
            T hz2inv = one/(dz*dz);
            T dxinv = one/dx;
            T dzinv = one/dz;
            T a = hx2inv + hz2inv;
            T b = four*(dxinv*dt0dx + dzinv*dt0dz)
                - two*(tauX*hx2inv + tauZ*hz2inv);
            T c = (tauX*tauX)*hx2inv + (tauZ*tauZ)*hz2inv
                - four*((dt0dx*tauX)*dxinv + (dt0dz*tauZ)*dzinv)
                + four*(sourceSlowness*sourceSlowness - s0*s0); 
            T det = b*b - four*a*c;
//std::cout << "apoly: " << a << " " << b << " " << c << " " << t0 << std::endl;
            if (det >= 0)
            {
                T tau = (-b + sycl::sqrt(det))/(two*a);
                t0Update3 = t0 + tau; 
            }
            if (t0Update3 < t1 || t0Update3 < t3){t0Update3 = huge;}
//std::cout << t0Update3 << std::endl;
        }
        tUpd = sycl::fmin(t0Update1, t0Update3);
    }
    return tUpd;
}

/*
/// @brief The 2D Cartesian finite-difference operators for headwaves.
/// @param[in] dx   The grid spacing in x in m.
/// @param[in] dz   The grid spacing in z in m.
/// @param[in] s0   The slowness in the cell in s/m. 
/// @param[in] s1   The slowness in the cell and the cell to the right or
///                 left of s0 in s/m.
/// @param[in] s3   The slowness in the cell and the cell to the top or
///                 or bottom of s0 in s/m.
/// @param[in] t1   The travel time at the grid point to the left or
///                 right of the update node in seconds.
/// @param[in] t3   The travel time at the grid point to the top or
///                 bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
/// @note This is for the general case when the grid spacing in x does not
///       equal the grid spacing in z.
#pragma omp declare simd uniform(dx, dz)
template<typename T>
[[nodiscard]]
T cartesianFiniteDifferenceTwoPoint(const T dx, const T dz,
                                    const T s0, const T s1, const T s3,
                                    const T t1, const T t3)
{
    T t0x = t1 + dx*sycl::fmin(s0, s3);
    T t0z = t3 + dz*sycl::fmin(s0, s1);
    return sycl::fmin(t0x, t0z);
}

/// @brief The 2D Cartesian finite-difference operators for headwaves.
/// @param[in] h    The grid spacing in x and z in m.
/// @param[in] s0   The slowness in the cell in s/m.
/// @param[in] s1   The slowness in the cell and the cell to the right or
///                 left of s0 in s/m.
/// @param[in] s3   The slowness in the cell and the cell to the top or
///                 or bottom of s0 in s/m.
/// @param[in] t1   The travel time at the grid point to the left or
///                 right of the update node in seconds.
/// @param[in] t3   The travel time at the grid point to the top or
///                 bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
#pragma omp declare simd uniform(h)
template<typename T>
[[nodiscard]]
T cartesianFiniteDifferenceTwoPoint(const T h,
                                    const T s0, const T s1, const T s3,
                                    const T t1, const T t3)
{
    T t0x = t1 + h*sycl::fmin(s0, s3);
    T t0z = t3 + h*sycl::fmin(s0, s1);
    return sycl::fmin(t0x, t0z);
}

/// @brief The 2D Cartesian 3-point finite-difference operators.
/// @param[in] huge      A default travel time (seconds) should the illumination
///                      conditions not be satisfied.
/// @param[in] dx        The grid spacing in x in m.
/// @param[in] dz        The grid spacing in z in m.
/// @param[in] dx_dz     dx/dz
/// @param[in] dz_dx     dz/dx
/// @param[in] cosTheta  dx/(dx**2 + dz**2)
/// @param[in] sinTheta  dz/(dx**2 + dz**2)
/// @param[in] s0        The slowness in the cell in s/m. 
/// @param[in] t1        The travel time at the grid point to the left or
///                      right of the update node in seconds.
/// @param[in] t2        The travel time at the grid point in opposing
///                      corner of the current cell in seconds.
/// @param[in] t3        The travel time at the grid point to the top or
///                      bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
/// @note This is for the general case when the grid spacing in x does not
///       equal the grid spacing in z.
#pragma omp declare simd uniform(huge, dx, dz, dx2i, dz2i, cosTheta, sinTheta)
template<typename T>
[[nodiscard]]
T cartesianFiniteDifferenceThreePoint(const T huge,
                                      const T dx, const T dz,
                                      const T dx_dz, const T dz_dx,
                                      const T cosTheta, //dx/sqrt(dx*dx + dz*dz)
                                      const T sinTheta, //dz/sqrt(dx*dx + dz*dz)
                                      const T s0,
                                      const T t1, const T t2, const T t3)
{
    T s02 = s0*s0;
    T hzs0 = dz*s0;
    T dTz = t1 - t2;
    T sz = s0*sinTheta;
    T t0a = huge;
    T deta = hzs0*hzs0 - dTz*dTz;
    if (dTz >= 0 && deta >= 0 && dTz <= dz*sz){t0a = t1 + dx_dz*sycl::sqrt(deta);}

    T hxs0 = dx*s0;
    T dTx = t3 - t2;
    T sx = s0*cosTheta;
    T t0b = huge;
    T detb = hxs0*hxs0 - dTx*dTx; 
    if (dTx >= 0 && detb >= 0 && dTx <= dz*sx){t0b = t1 + dz_dx*sycl::sqrt(detb);}
    return sycl::fmin(t0a, t0b);
}

/// @brief The 2D Cartesian 3-point finite-difference operators.
/// @param[in] huge   A default travel time (seconds) should the illumination
///                   conditions not be satisfied.
/// @param[in] h      The grid spacing in x and z in m.
/// @param[in] s0     The slowness in the cell in s/m.
/// @param[in] t1     The travel time at the grid point to the left or
///                   right of the update node in seconds.
/// @param[in] t2     The travel time at the grid point in opposing
///                   corner of the current cell in seconds.
/// @param[in] t3     The travel time at the grid point to the top or
///                   bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
#pragma omp declare simd uniform(huge, h)
template<class T>
[[nodiscard]]
T cartesianFiniteDifferenceThreePoint(const T huge, const T h,
                                      const T s0,
                                      const T t1, const T t2, const T t3)
{
    T hs = h*s0;

    T dTz = t1 - t2; 
    T t0a = huge;
    T deta = hs*hs - dTz*dTz;
    if (dTz >= 0 && deta >= 0 && dTz <= M_SQRT1_2*hs)
    {
        t0a = t1 + sycl::sqrt(deta);
    }

    T dTx = t3 - t2;
    T t0b = huge;
    T detb = hs*hs - dTx*dTx;
    if (dTx >= 0 && detb >= 0 && dTx <= M_SQRT1_2*hs)
    {
        t0b = t2 + sycl::sqrt(detb);
    }
    return sycl::fmin(t0a, t0b);

}
/// @brief The 2D Cartesian four-point finite-difference operators.
/// @param[in] huge           The default travel time should the illumination
///                           conditions not be satisfied.
/// @param[in] dx             The grid spacing in x in m.
/// @param[in] dz             The grid spacing in z in m.
/// @param[in] dx2            dx^2 in m^2.
/// @param[in] dz2            dz^2 in m^2.
/// @param[in] dx2_p_dz2_inv  1/(dx^2 + dz^2) in 1/m^2.
/// @param[in] s0             The slowness in the cell in s/m.
/// @param[in] t1             The travel time at the grid point to the left or
///                           right of the update node in seconds.
/// @param[in] t2             The travel time at the grid point to the
///                           far corner of the update node in seconds.
/// @param[in] t3             The travel time at the grid point to the top or
///                           bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
/// @note This is for the general case when the grid spacing in x does not
///       equal the grid spacing in z.
#pragma omp declare simd uniform(huge, dx, dz, dx2, dz2, dx2_p_dz2_inv)
template<typename T>
[[nodiscard]]
T carestianFiniteDifferenceFourPoint(const T huge,
                                     const T dx, const T dz,
                                     const T dx2, const T dz2,
                                     const T dx2_p_dz2_inv,
                                     const T s0,
                                     const T t1, const T t2, const T t3)
{
    T dTx = t3 - t2;
    T dTz = t1 - t2;
    T Tx = t1 - dTx; //t1 - (t3 - t2) = t1 - t3 + t2
    T Tz = t3 - dTz; //t3 - (t1 - t2) = t3 - t1 + t2
    T T1mT3 = t1 - t3;
    T det = (dx2 + dz2)*(s0*s0) - T1mT3*T1mT3;
    // Require the wave is propagating in the right direction.
    // Additionally, only real roots are found when det >= 0.
    T t0 = huge;
    if (dTx >= 0 && dTz >= 0 && det >= 0 && dTx <= dx*s0 && dTz <= dz*s0)
    {
        t0 = dx2_p_dz2_inv*((Tz*dx2 + Tx*dx2) + 2*(dx*dz)*sycl::sqrt(det));
    }
    if (t0 < 0){t0 = huge;}
    return t0;
}

/// @brief The 2D Cartesian four-point finite-difference operators.
/// @param[in] huge  The default travel time should the illumination conditions
///                  not be satisfied.
/// @param[in] h     The grid spacing in x and z in m.
/// @param[in] s0    The slowness in the cell in s/m.
/// @param[in] t1    The travel time at the grid point to the left or
///                  right of the update node in seconds.
/// @param[in] t2    The travel time at the grid point to the
///                  far corner of the update node in seconds.
/// @param[in] t3    The travel time at the grid point to the top or
///                  bottom of the update node in seconds.
/// @result The candidate update travel time in seconds.
#pragma omp declare simd uniform(huge, h)
template<typename T>
[[nodiscard]]
T carestianFiniteDifferenceFourPoint(const T huge, const T h,
                                     const T s0,
                                     const T t1, const T t2, const T t3)
{
    T hs = h*s0;
    T dTx = t3 - t2;
    T dTz = t1 - t2;
    T T1mT3 = t1 - t3;
    T det = 2*(hs*hs) - T1mT3*T1mT3; 
    T t0 = huge;
    if (dTx >= 0 && dTz >= 0 && det >= 0 && dTx <= hs && dTz <= hs)
    {
        t0 = t2 + sycl::sqrt(det);
    }
    return t0;
}


/// @brief The Cartesian finite difference.
/// @param[in] dx   The grid spacing in x in m.
/// @param[in] dz   The grid spacing in z in m.
#pragma omp declare simd uniform(dx, dz)
template<typename T>
T cartesianFiniteDifference(const T huge,
                            const T dx, const T dz,
                            const T dx2i, const T dz2i,
                            const T cosTheta, const T sinTheta,
                            const T *s, const T *t)
{
    /// 1D operators for headwaves
    auto t1 = cartesianFiniteDifferenceTwoPoint(dx, dz,
                                                s[0], s[1], s[3],
                                                t[1], t[3]);
    // 2D three-point operators
    auto t2 = cartesianFiniteDifferenceThreePoint(huge,
                                                  dx, dz,
                                                  dx2i, dz2i,
                                                  cosTheta, sinTheta,
                                                  s[0],
                                                  t[1], t[2], t[3]);
    // 2D four-point operators
    auto t3 = carestianFiniteDifferenceFourPoint(huge,
                                                 dx, dz,
                                                 dx2i, dz2i,
                                                 s[0],
                                                 t[1], t[2], t[3]); 
    /// 
    return sycl::fmin(sycl::fmin(t1, t2), t3);
}
*/

/*
/// @brief Converts a grid index (ix,iz) to an index.
/// @param[in] nx   The number of grid point in x.
/// @param[in] ix   The ix'th grid point.
/// @param[in] iz   The iz'th grid point.
/// @result The index in a row-major [nz x nx] matrix corresponding to the
///         given (ix, iz)'th grid point.
#pragma omp declare(simd) uniform(nx)
int gridToIndex(const int nx, const int ix, const int iz)
{
    return iz*nx + ix;
}
*/

/// @brief Gets the sweep finite difference signs and shifts.
/// @param[out] ixShift  The shift in the x grid point when computing the
///                      analytic travel times and derivatives.
/// @param[out] izShift  The shift in the z grid point when computing the
///                      analytic travel times and derivatives.
/// @param[out] signX    The sign on the derivative in x.  This makes the
///                      derivative sign consistent with the sweep direction.
/// @param[out] signZ    The sign on the derivative in z.  This makes the
///                      derivative sign consistent with the sweep direction.
template<EikonalXX::SweepNumber2D E>
void getSweepFiniteDifferenceSigns(int *ixShift, int *izShift,
                                   int *signX, int *signZ)
{
    *ixShift =-1; 
    *izShift =-1; 
    *signX = 1; 
    *signZ = 1; 
    if constexpr (E == SweepNumber2D::SWEEP2)
    {
        *signX =-1; 
        *ixShift = 1; 
    }
    else if constexpr (E == SweepNumber2D::SWEEP3)
    {
        *signZ =-1; 
        *izShift = 1; 
    }
    else if constexpr (E == SweepNumber2D::SWEEP4)
    {
        *signX =-1; 
        *signZ =-1; 
        *ixShift = 1; 
        *izShift = 1; 
    }
}

/// @brief Gets the loop limits for the fast sweeping method i.e.,:
///        for (int iz = iz0; iz != iz1; iz = iz + izDir)
///            for (int ix = ix0; ix != ix1; ix = ix + ixDir)
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridZ  The number of grid points in z.
/// @param[out] ix0    The first x grid index in the loop.
/// @param[out] ix1    The stopping x grid index (exclusive).
/// @param[out] iz0    The first z grid index in the loop.
/// @param[out] iz1    The stopping z grid index (exclusive).
/// @param[out] ixDir  The x loop variable update increment (+1 or -1).
/// @param[out] izDir  The z loop variable update increment (+1 or -1).
template<EikonalXX::SweepNumber2D E>
void getLoopLimits(const int nGridX, const int nGridZ,
                   int *ix0, int *iz0,
                   int *ix1, int *iz1,
                   int *ixDir, int *izDir)
{
    if constexpr (E == SweepNumber2D::SWEEP1)
    {
        *ix0 = 1;
        *ix1 = nGridX;
        *iz0 = 1;
        *iz1 = nGridZ;
        *ixDir = 1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber2D::SWEEP2)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iz0 = 1;
        *iz1 = nGridZ;
        *ixDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber2D::SWEEP3)
    {
        *ix0 = 1;
        *ix1 = nGridX;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir = 1;
        *izDir =-1;
    }
    else //if (sweep == SweepNumber2D::SWEEP4)
    {
        *ix0 = nGridX - 2;
        *ix1 =-1;
        *iz0 = nGridZ - 2;
        *iz1 =-1;
        *ixDir =-1;
        *izDir =-1;
    }
}

/// @brief Gets the loop limits for the initialization of the
///        fast sweeping method.  Effectively, during initialization
///        we only care about propagating energy from the surface
///        and leave it to the refinement loops to clean up the
///        solution 
/// @param[in] iSrcX   The source grid point in x.
/// @param[in] iSrcZ   The source grid point in z.
/// @param[in] nGridX  The number of grid points in x.
/// @param[in] nGridZ  The number of grid points in z.
/// @param[out] ix0    The first x grid index in the loop.
/// @param[out] ix1    The stopping x grid index (exclusive).
/// @param[out] iz0    The first z grid index in the loop.
/// @param[out] iz1    The stopping z grid index (exclusive).
/// @param[out] ixDir  The loop variable update increment (+1 or -1).
/// @param[out] izDir  The loop variable update increment (+1 or -1).
template<EikonalXX::SweepNumber2D E>
void getLoopLimits(const int iSrcX, const int iSrcZ,
                   const int nGridX, const int nGridZ,
                   int *ix0, int *iz0,
                   int *ix1, int *iz1,
                   int *ixDir, int *izDir)
{
    if constexpr (E == SweepNumber2D::SWEEP1)
    {
        *ix0 = std::max(iSrcX, 1);
        *ix1 = nGridX;
        *iz0 = std::max(iSrcZ, 1);
        *iz1 = nGridZ;
        *ixDir = 1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber2D::SWEEP2)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iz0 = std::max(iSrcZ, 1);
        *iz1 = nGridZ;
        *ixDir =-1;
        *izDir = 1;
    }
    else if constexpr (E == SweepNumber2D::SWEEP3)
    {
        *ix0 = std::max(iSrcX, 1);
        *ix1 = nGridX;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
        *iz1 =-1;
        *ixDir = 1;
        *izDir =-1;
    }
    else //if (sweep == SweepNumber2D::SWEEP4)
    {
        *ix0 = std::min(iSrcX + 1, nGridX - 2);
        *ix1 =-1;
        *iz0 = std::min(iSrcZ + 1, nGridZ - 2);
        *iz1 =-1;
        *ixDir =-1;
        *izDir =-1;
    }
}


/// @brief Converts a (sweep, level, index) to a grid point (ix, iz).
/// @param[in] sweep  The sweep number.
/// @param[in] level  The level in the level set method.
/// @param[in] indx   The index in the level.
/// @param[in] nx     The number of grid points in x.
/// @param[in] nz     The number of grid points in z.
/// @param[out] ix    The corresponding grid index in x.
/// @param[out] iz    The corresponding grid index in z.
#pragma omp declare simd uniform(sweep, level, nx, nz)
void sweepLevelIndexToGrid(const SweepNumber2D sweep, //const int sweep,
                           const int level,
                           const int indx,
                           const int nx,
                           const int nz,
                           int *ix, int *iz)
{
    // +x and +z 
    *iz = indx;
    *ix = level - *iz;
    // -x and +z
    if (sweep == SweepNumber2D::SWEEP2)
    {
        *ix = nx - 1 - *ix;
    }
    else if (sweep == SweepNumber2D::SWEEP3) // +x and -z
    {
        *iz = nz - 1 - *iz;
    }
    else if (sweep == SweepNumber2D::SWEEP4) // -x and -z
    {
        *iz = nz - 1 - *iz;
        *ix = nx - 1 - *ix;
    }
}

/// @brief Converts a grid point (ix, iz) in a sweep to the corresponding
///        level and index.
/// @param[in] sweep   The sweep number.
/// @param[in] ix      The grid index in x.
/// @param[in] iz      The grid index in z.
/// @param[in] nx      The number of grid points in x.
/// @param[in] nz      The number of grid points in z.
/// @param[out] level  The corresponding level.
/// @param[out] indx   The corresponding indx in the level.  Realize, this does
///                    not necessarily begin at 0.  To obtain this information
///                    use \c getLevelStartStopIndices() and subtract by i0.
/// @note This is the inverse operation of \c sweepLevelIndexToGrid().
#pragma omp declare simd uniform(sweep, nx ,nz)
void gridSweepToLevelIndex(const SweepNumber2D sweep, //const int sweep,
                           const int ix,
                           const int iz,
                           const int nx,
                           const int nz,
                           int *level, int *indx)
{
    auto jx = ix;
    auto jz = iz;
    if (sweep == SweepNumber2D::SWEEP2)
    {
        jx = nx - 1 - ix;
    }
    else if (sweep == SweepNumber2D::SWEEP3)
    {
        jz = nz - 1 - iz;
    }
    else if (sweep == SweepNumber2D::SWEEP4)
    {
        jz = nz - 1 - iz;
        jx = nx - 1 - ix;
    } 
    *indx = jz;
    *level = jx + jz;
}

/// @brief Converts a (sweep, level, index) to a travel time field index.
/// @param[in] sweep  The sweep number.
/// @param[in] level  The level in the level set method.
/// @param[in] indx   The index in the level.
/// @param[in] nx     The number of grid points in x.
/// @param[in] nz     The number of grid points in z.
/// @result The index in the travel time field corresponding to the sweep, level,
///                   and index in the level.
#pragma omp declare simd uniform(sweep, level, nx, nz)
[[maybe_unused]]
int sweepLevelIndexToIndex(const SweepNumber2D sweep, //const int sweep,
                           const int level,
                           const int indx,
                           const int nx,
                           const int nz)
{
    int ix, iz;
    sweepLevelIndexToGrid(sweep, level, indx, nx, nz, &ix, &iz);
    return gridToIndex(nx, ix, iz);
}

/// @brief Converts a (sweep, level, index) to the grid points for finite
///        differencing a travel time.
/// @param[in] sweep  The sweep number.
/// @param[in] level  The level in the level set method.
/// @param[in] indx   The index in the level.
/// @param[in] nx     The number of grid points in x.
/// @param[in] nz     The number of grid points in z.
/// @param[out] it0   The corresponding grid index at the update node.
/// @param[out] it1   The corresponding grid index to the left/right of the
///                   update node.
/// @param[out[ it2   The corresponding grid index to the across from the
///                   update node. 
/// @param[out] it3   The corresponding grid index above/below the update node.
#pragma omp declare simd uniform(sweep, level, nx, nz)
[[maybe_unused]]
void sweepLevelIndexToTravelTimeIndices(
    const SweepNumber2D sweep, const int level, const int indx,
    const int nx, const int nz,
    int *it0, int *it1, int *it2, int *it3)
{
    int ix, iz;
    sweepLevelIndexToGrid(sweep, level, indx, nx, nz, &ix, &iz);
    *it0 = gridToIndex(nx, ix, iz);
    if (sweep == SweepNumber2D::SWEEP1)
    {
        *it1 = gridToIndex(nx, std::max(0, ix-1), iz); 
        *it2 = gridToIndex(nx, std::max(0, ix-1), std::max(0, iz-1));
        *it3 = gridToIndex(nx, ix,                std::max(0, iz-1)); 
    }
    else if (sweep == SweepNumber2D::SWEEP2)
    {
        *it1 = gridToIndex(nx, std::min(nx-1, ix+1), iz);
        *it2 = gridToIndex(nx, std::min(nx-1, ix+1), std::max(0, iz-1));
        *it3 = gridToIndex(nx, ix,                   std::max(0, iz-1));
    } 
    else if (sweep == SweepNumber2D::SWEEP3)
    {
        *it1 = gridToIndex(nx, std::max(0, ix-1), iz);
        *it2 = gridToIndex(nx, std::max(0, ix-1), std::min(nz-1, iz+1));
        *it3 = gridToIndex(nx, ix,                std::min(nz-1, iz+1));
    }
    else // sweep == SweepNumber2D::SWEEP4 
    {
        *it1 = gridToIndex(nx, std::min(nx-1, ix+1), iz);
        *it2 = gridToIndex(nx, std::min(nx-1, ix+1), std::min(nz-1, iz+1));
        *it3 = gridToIndex(nx, ix,                   std::min(nz-1, iz+1));
    }
}

/// @brief For each level this returns the start/stop indices in the foor loop
///        (e.g.,  for (int i=i0; i<i1; ++i){}).
/// @param[in] nx     The number of x grid points.
/// @param[in] nz     The number of z grid points.
/// @param[in] level  The current level.  This is C indexed.
/// @param[out] i0    The level start index (inclusive)
/// @param[out] i1    The level stop idnex (exclusive).
void getLevelStartStopIndices(const int nx, const int nz, const int level,
                              int *i0, int *i1) noexcept
{
     *i0 = sycl::max(0, level + 1 - nx);
     *i1 = sycl::min(nz, level + 1);
}

/// @brief From the number of grid points this returns the number of levels.
/// @param[in] nx  The number of grid points in x.
/// @param[in] nz  The number of grid points in z.
/// @result The number of levels in the 2D solver.
int computeNumberOfLevels(const int nx, const int nz) noexcept
{
    return nx + nz - 1;
}

/// @brief During the solution phase we need to point to the start of a 
///        level-set method specific start index in the travel time field.
///        The vector generated by this function makes that straightforward
///        to do.
/// @param[in] nx  The number of grid points in x.
/// @param[in] nz  The number of grid points in z.
/// @result The nodes in the travel time field for the level set method
///         at the level'th level begin at levelOffset[level].
std::vector<int> makeLevelOffset(const int nx, const int nz)
{
    auto nLevels = computeNumberOfLevels(nx, nz);
    std::vector<int> levelOffset(nLevels + 1, -1);
    int i0, i1;
    levelOffset[0] = 0;
    for (int level=0; level<nLevels; ++level)
    {
        getLevelStartStopIndices(nx, nz, level, &i0, &i1);
        levelOffset[level+1] = levelOffset[level] + (i1 - i0);
    }
    return levelOffset;
}

/// @brief Converts the grid point to the surrounding slowness cells.
/// @param[in] sweep    The sweep number.
/// @param[in] ix       The ix'th grid point.
/// @param[in] iz       The iz'th grid point.
/// @param[in] nCellX   The number of cells in x.
/// @param[in] nCellZ   The number of cells in z.
/// @param[out] iCell0  The index of the slowness field corresponding to 
///                     the slowness of the home cell.
/// @param[out] iCell1  The index of the slowness field corresponding to
///                     the slowness to the left or right of the home cell.
/// @param[out] iCell3  The index of the slowness field corresponding to 
///                     the slowness above or below the home cell. 
void gridToSurroundingSlowness(const SweepNumber2D sweep, //int sweep,
                               const int ix, const int iz,
                               const int nCellX, const int nCellZ,
                               int *iCell0, int *iCell1, int *iCell3)
{
    int iCell0X, iCell0Z, iCell1X, iCell1Z,
        iCell2X, iCell2Z, iCell3X, iCell3Z = 0;
    if (sweep == SweepNumber2D::SWEEP1)
    {
        iCell0X = sycl::max(0, ix - 1);
        iCell1X = sycl::min(nCellX - 1, iCell0X + 1);
        if (ix == 0){iCell1X = 0;}
        iCell2X = iCell1X;
        iCell3X = iCell0X;

        iCell0Z = sycl::max(0, iz - 1);
        iCell2Z = sycl::min(nCellZ - 1, iCell0Z + 1);
        if (iz == 0){iCell2Z = 0;}
        iCell1Z = iCell0Z;
        iCell3Z = iCell2Z;
    }
    else if (sweep == SweepNumber2D::SWEEP2)
    {
        iCell0X = ix;
        iCell1X = sycl::max(0, iCell0X - 1);
        if (ix == nCellX)
        {
            iCell0X = ix - 1;
            iCell1X = ix - 1;
        }
        iCell2X = iCell1X;
        iCell3X = iCell0X;

        iCell0Z = sycl::max(0, iz - 1);
        iCell2Z = sycl::min(nCellZ - 1, iCell0Z + 1);
        if (iz == 0){iCell2Z = 0;}
        iCell1Z = iCell0Z;
        iCell3Z = iCell2Z;
    }
    else if (sweep == SweepNumber2D::SWEEP3)
    {
        iCell0X = sycl::max(0, ix - 1);
        iCell1X = sycl::min(nCellX - 1, iCell0X + 1);
        if (ix == 0){iCell1X = 0;}
        iCell2X = iCell1X;
        iCell3X = iCell0X;

        iCell0Z = iz;
        iCell2Z = sycl::max(0, iCell0Z - 1);
        if (iz == nCellZ)
        {
            iCell0Z = iz - 1;
            iCell2Z = iz - 1;
        }
        iCell1Z = iCell0Z;
        iCell3Z = iCell2Z;
    }
    else //if (sweep == SweepNumber2D::SWEEP4)
    {
        iCell0X = ix;
        iCell1X = sycl::max(0, iCell0X - 1);
        if (ix == nCellX)
        {
            iCell0X = ix - 1;
            iCell1X = ix - 1;
        }
        iCell2X = iCell1X;
        iCell3X = iCell0X;

        iCell0Z = iz;
        iCell2Z = sycl::max(0, iCell0Z - 1);
        if (iz == nCellZ)
        {
            iCell0Z = iz - 1;
            iCell2Z = iz - 1;
        }
        iCell1Z = iCell0Z;
        iCell3Z = iCell2Z;
    } 
    *iCell0 = gridToIndex(nCellX, iCell0X, iCell0Z);
    *iCell1 = gridToIndex(nCellX, iCell1X, iCell1Z);
    //*iCell2 = gridToIndex(nCellX, iCell2X, iCell2Z);
    *iCell3 = gridToIndex(nCellX, iCell3X, iCell3Z);
}

/// @brief Converts a grid point to surrounding travel times indices.
/// @param[in] ix      The grid point in x.
/// @param[in] iz      The grid point in z.
/// @param[in] nGridX  The number of grid points in x. 
/// @param[in] it0     The travel time field index of the update node.
/// @param[in] it1     The travel time field index to the left or right of th
///                    update node.
/// @param[in] it2     The travel time field index of diagonal to the
///                    update node.
/// @param[in] it3     The travel time field index above or below the
///                    update node.
void gridToSurroundingTravelTimes(const SweepNumber2D sweep, //int sweep,
                                  const int ix, const int iz,
                                  const int nGridX,
                                  int *it0, int *it1, int *it2, int *it3)
{
    if (sweep == SweepNumber2D::SWEEP1)
    {
        *it0 = gridToIndex(nGridX, ix,     iz);
        *it1 = gridToIndex(nGridX, ix - 1, iz);
        *it2 = gridToIndex(nGridX, ix - 1, iz - 1);
        *it3 = gridToIndex(nGridX, ix,     iz - 1);
    }
    else if (sweep == SweepNumber2D::SWEEP2)
    {
        *it0 = gridToIndex(nGridX, ix,     iz);
        *it1 = gridToIndex(nGridX, ix + 1, iz);
        *it2 = gridToIndex(nGridX, ix + 1, iz - 1);
        *it3 = gridToIndex(nGridX, ix,     iz - 1);
    }
    else if (sweep == SweepNumber2D::SWEEP3)
    {
        *it0 = gridToIndex(nGridX, ix,     iz);
        *it1 = gridToIndex(nGridX, ix - 1, iz);
        *it2 = gridToIndex(nGridX, ix - 1, iz + 1);
        *it3 = gridToIndex(nGridX, ix,     iz + 1);
    }
    else // (sweep == SweepNumber2D::SWEEP4) 
    {
        *it0 = gridToIndex(nGridX, ix,     iz);
        *it1 = gridToIndex(nGridX, ix + 1, iz);
        *it2 = gridToIndex(nGridX, ix + 1, iz + 1);
        *it3 = gridToIndex(nGridX, ix,     iz + 1);
    }
}
 

/// @brief Converts the slowness model to a slowness model for the sweep.
/// @param[in] sweep        The current sweep number.
/// @param[in] nLevels      The number of levels in the solver.
/// @param[in] nGridX       The number of x grid points.
/// @param[in] nGridZ       The number of z grid points.
/// @param[in] nCellX       The number of cells in x (= nGridX - 1).
/// @param[in] nCellZ       The number of cells in z (= nGridZ - 1).
/// @param[in] slow         The slowness field.  This is an [nCellZ x nCellX]
///                         matrix in row major format.
/// @param[out] sweepSlowness  The slownesses surrounding each grid point in
///                            the level-set solver. 
template<class T>
void slownessToSweepSlowness(const SweepNumber2D sweep, //const int sweep,
                             const size_t nLevels,
                             const int nGridX, const int nGridZ,
                             const int nCellX, const int nCellZ,
                             const T *slow,
                             SweepSlowness2D<T> *sweepSlowness)
{
    //auto nCell = static_cast<size_t> (nCellX)*static_cast<size_t> (nCellZ); 
    //auto nGrid = static_cast<size_t> (nGridX)*static_cast<size_t> (nGridZ);
    tbb::parallel_for(tbb::blocked_range<int> (0, nLevels),
                      [=](const tbb::blocked_range<int> &localLevel)
    {
        int iDst, i0, i1, ix, iz;
        int iCell0, iCell1, iCell3 = 0;
        //int iCell0X, iCell1X, iCell2X, iCell3X = 0;
        //int iCell0Z, iCell1Z, iCell2Z, iCell3Z = 0;
        for (int level=localLevel.begin(); level != localLevel.end(); ++level)
        {
            auto s0 = sweepSlowness[level].s0;
            auto s1 = sweepSlowness[level].s1;
            auto s3 = sweepSlowness[level].s3;
            getLevelStartStopIndices(nGridX, nGridZ, level, &i0, &i1);
            for (int indx=i0; indx<i1; ++indx)
            {
                // Get (ix, iz) grid node
                sweepLevelIndexToGrid(sweep, level, indx, nGridX, nGridZ,
                                      &ix, &iz);
                // Get the slowness neighbors
                gridToSurroundingSlowness(sweep,
                                          ix, iz,
                                          nCellX, nCellZ,
                                          &iCell0, &iCell1, &iCell3);
                iDst = indx - i0;
                s0[iDst] = slow[iCell0];
                s1[iDst] = slow[iCell1];
                s3[iDst] = slow[iCell3];
            } // Loop on nodes in sweep
        } // Loop on levels
    }); // End lambda
}

/// @brief Sets the preliminary nodes available for updating based on a sweep.
/// @param[in] sweep         The sweep number.
/// @param[in] nLevels       The number of levels in the solver.
/// @param[in] nGridX        The number of x grid points.
/// @param[in] nGridZ        The number of z grid points.
/// @param[in] levelOffset   A map from the level to the offset in the
///                          travel time field.  This is an array with
///                          dimension [nLevels + 1].
/// @param[out] lUpdateNode  UPDATE_NODE indicates that this node in the level
///                          can update the travel time field.
///                          BOUNDARY_NODE indicates that this node cannot
///                          be updated.
void setPreliminaryUpdateNodes(const SweepNumber2D sweep, //const int sweep,
                               const size_t nLevels,
                               const int nGridX, const int nGridZ,
                               const int *levelOffset,
                               int8_t *lUpdateNode)
{
    tbb::parallel_for(tbb::blocked_range<int> (0, nLevels),
                      [=](const tbb::blocked_range<int> &localLevel)
    {
        for (int level=localLevel.begin(); level != localLevel.end(); ++level)
        {
            int i0, i1, ix, iz;
            auto offset = levelOffset[level];
            getLevelStartStopIndices(nGridX, nGridZ, level, &i0, &i1);
            for (int indx=i0; indx<i1; ++indx)
            {
                // Get (ix, iz) grid node
                sweepLevelIndexToGrid(sweep, level, indx, nGridX, nGridZ,
                                      &ix, &iz);                
                int8_t lUpdate = UPDATE_NODE;
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    if (ix == 0 || iz == 0){lUpdate = BOUNDARY_NODE;}
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    if (ix == nGridX - 1 || iz == 0){lUpdate = BOUNDARY_NODE;}
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    if (ix == 0 || iz == nGridZ - 1){lUpdate = BOUNDARY_NODE;}
                } 
                else //if (sweep == SweepNumber2D::SWEEP4)
                {
                    if (ix == nGridX - 1 || iz == nGridZ - 1)
                    {
                        lUpdate = BOUNDARY_NODE;
                    }
                }
                auto iDst = offset + (indx - i0);
                lUpdateNode[iDst] = lUpdate;
           } // Loop on nodes in level
        } // Loop on processes levels
    }); // End TBB parallel
}

/// @brief Sets the preliminary nodes available for updating based on a sweep
///        in the fast sweeping method.
/// @param[in] sweep         The sweep number.
/// @param[in] nGridX        The number of grid points in x.
/// @param[in] nGridZ        The number of grid points in z.
/// @param[out] lUpdateNode  UPDATE_NODE indicates this grid point can be
///                          updated.
///                          BOUNDARY_NODE indicates this grid point cannot
///                          be updated.
void setPreliminaryUpdateNodes(const SweepNumber2D sweep,
                               const int nGridX, const int nGridZ,
                               int8_t *lUpdateNode)
{
    auto nGrid = nGridX*nGridZ;
    std::fill(lUpdateNode, lUpdateNode + nGrid, UPDATE_NODE);
    int ix, iz;
    if (sweep == SweepNumber2D::SWEEP1)
    {
        iz = 0;
        #pragma omp simd
        for (ix = 0; ix < nGridX; ++ix)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
        ix = 0;
        #pragma omp simd
        for (iz = 0; iz < nGridZ; ++iz)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
    }
    else if (sweep == SweepNumber2D::SWEEP2)
    {
        iz = 0;
        #pragma omp simd
        for (ix = 0; ix < nGridX; ++ix)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
        ix = nGridX - 2;
        #pragma omp simd
        for (iz = 0; iz < nGridZ; ++iz)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
    }
    else if (sweep == SweepNumber2D::SWEEP3)
    {
        iz = nGridZ - 2;
        #pragma omp simd
        for (ix = 0; ix < nGridX; ++ix)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
        ix = 0;
        #pragma omp simd
        for (iz = 0; iz < nGridZ; ++iz)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
    }
    else //if (sweep == SweepNumber2D::SWEEP4)
    {
        iz = nGridZ - 2;
        #pragma omp simd
        for (ix = 0; ix < nGridX; ++ix)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
        ix = nGridX - 2;
        #pragma omp simd
        for (iz = 0; iz < nGridZ; ++iz)
        {
            lUpdateNode[gridToIndex(nGridX, ix, iz)] = BOUNDARY_NODE;
        }
    }
}
/*
    sycl::queue q{sycl::cpu_selector{}};
    {
    sycl::buffer<T> sweepSlowBuffer(sycl::range{4*nGrid});
    sweepSlowBuffer.set_final_data(sweepSlowness);
    q.submit([&](sycl::handler &h)
    {
        sycl::accessor slowAcc(slowBuffer, h, sycl::read_only);
        sycl::accessor sweepSlowAcc(sweepSlowBuffer, h, sycl::write_only);
        h.parallel_for(nLevels, [=](sycl::id<1> level)
        {
            int i0, i1, ix, iz;
            int iCell0X, iCell1X, iCell2X, iCell3X = 0;
            int iCell0Z, iCell1Z, iCell2Z, iCell3Z = 0;
            auto offset = levelOffset[level];
            getLevelStartStopIndices(nGridX, nGridZ, level, &i0, &i1);
            for (int indx=i0; indx<i1; ++indx)
            {
                // Get (ix, iz) grid node
                sweepLevelIndexToGrid(sweep, level, indx, nGridX, nGridZ,
                                      &ix, &iz);
                // Get the neighbors surrounding the grid point
                if (sweep == 0)
                {
                    iCell0X = sycl::max(0, ix - 1);
                    iCell1X = sycl::min(nCellX - 1, iCell0X + 1);
                    if (ix == 0){iCell1X = 0;}
                    iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = sycl::max(0, iz - 1);
                    iCell2Z = sycl::min(nCellZ - 1, iCell0Z + 1);
                    if (iz == 0){iCell2Z = 0;}
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }
                else if (sweep == 1)
                {
                    iCell0X = ix;
                    iCell1X = sycl::max(0, iCell0X - 1);
                    if (ix == nCellX)
                    {
                        iCell0X = ix - 1;
                        iCell1X = ix - 1;
                    }
                    iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = sycl::max(0, iz - 1);
                    iCell2Z = sycl::min(nCellZ - 1, iCell0Z + 1);
                    if (iz == 0){iCell2Z = 0;}
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z; 
                }
                else if (sweep == 2)
                {
                    iCell0X = sycl::max(0, ix - 1);
                    iCell1X = sycl::min(nCellX - 1, iCell0X + 1);
                    if (ix == 0){iCell1X = 0;}
                    iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = iz;
                    iCell2Z = sycl::max(0, iCell0Z - 1);
                    if (iz == nCellZ)
                    {
                        iCell0Z = iz - 1;
                        iCell2Z = iz - 1;
                    }
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }
                else // sweep == 3
                {
                    iCell0X = ix;
                    iCell1X = sycl::max(0, iCell0X - 1);
                    if (ix == nCellX)
                    {
                        iCell0X = ix - 1;
                        iCell1X = ix - 1;
                    }
                    iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = iz;
                    iCell2Z = sycl::max(0, iCell0Z - 1);
                    if (iz == nCellZ)
                    {
                        iCell0Z = iz - 1;
                        iCell2Z = iz - 1;
                    }
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }
                auto iCell0 = gridToIndex(nCellX, iCell0X, iCell0Z);
                auto iCell1 = gridToIndex(nCellX, iCell1X, iCell1Z);
                auto iCell2 = gridToIndex(nCellX, iCell2X, iCell2Z);
                auto iCell3 = gridToIndex(nCellX, iCell3X, iCell3Z);
                auto iDst0 = 4*(offset + (indx - i0));
                auto iDst1 = iDst0 + 1;
                auto iDst2 = iDst0 + 2;
                auto iDst3 = iDst0 + 3;
                sweepSlowAcc[iDst0] = slowAcc[iCell0];
                sweepSlowAcc[iDst1] = slowAcc[iCell1];
                sweepSlowAcc[iDst2] = slowAcc[iCell2];
                sweepSlowAcc[iDst3] = slowAcc[iCell3];
            } // Loop on points in level
        }); // End kernel
    }); // End q.submit
    // Copy result of buffer back
    //sycl::host_accessor sweepSlowAcc(sweepSlowBuffer);
    //for (size_t i=0; i<4*nCell; ++i)
    //{
    //    sweepSlowness[i] = sweepSlowAcc[i];
    //}
    } 
    q.wait();
    */

/// @brief A finite-difference's concept of `forward' and `backward' is
///        dependent on the sweep direction.  To access the appropriate
///        slowness cells and traveltime nodes for the sweep we use this
///        convenience function.
/// @param[in] sweep  The sweep number.  This must be 1, 2, 3, or 4. 
/// @param[out] shiftTx   The traveltime field offset in x for the local
///                       solver's finite difference.
/// @param[out] shiftTz   The traveltime field offset in z for the lcoal
///                       solver's finite difference.
/// @param[out] shiftSx   The slowness field offset in x for the local
///                       solver's finite difference.
/// @param[out] shiftSz   The slowness field offset in z for the local
///                       solver's finite difference.
/*
void getSweepFiniteDifferenceShifts(const int sweep,
                                    int *shiftTx, int *shiftTz,
                                    int *shiftSx, int *shiftSz) 
{
    if (sweep == 1) // +x and +z
    {
        *shiftTz = 1;
        *shiftTx = 1;
        *shiftSz = 1;
        *shiftSx = 1;
    }
    else if (sweep == 2) // -x and +z
    {
        *shiftTz = 1;
        *shiftTx =-1;
        *shiftSz = 1;
        *shiftSx = 0;
    }
    else if (sweep == 3) // +x and -z
    {
        *shiftTz =-1;
        *shiftTx = 1;
        *shiftSz = 0;
        *shiftSx = 1;
    }
    else if (sweep == 4) // -x and -z
    {
        *shiftTz =-1;
        *shiftTx =-1;
        *shiftSz = 0;
        *shiftSx = 0;
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}
*/

/*
/// @brief Performs the update.
template<class T>
class InitializedUpdate
{
public:
    InitializedUpdate(const T dx, const T dz,
                      const T dz2i_p_dx2i, const T dz2i_dx2i,
                      const T *localSlowness,
                      const T *travelTimes,
                      T *travelTimesUpdate) :
        mTravelTimes(travelTimes),
        mLocalSlowness(localSlowness),
        mTravelTimesUpdate(travelTimesUpdate),
        mDx(dx),
        mDz(dz)
    {
        auto dx2 = dx*dx;
        auto dz2 = dz*dz;
        mDx2i = one/dx2; 
        mDz2i = one/dz2;
        mDx2i_Dz2i = mDx2i*mDz2i;
        mDx2i_p_Dz2i = mDx2i + mDz2i;
        mDx2i_p_Dz2i_inv = one/mDx2i_p_Dz2i;
    }
    void operator()(const sycl::id<1> &i)
    {
        auto i4 = i*4; 
        // Extract slownesses
        auto s1 = mLocalSlowness[i4];
        auto s2 = mLocalSlowness[i4+1];
        auto s3 = mLocalSlowness[i4+2]; 
        // Extract travel times
        auto tv  = mTravelTimes[i4];
        auto te  = mTravelTimes[i4+1];
        auto tev = mTravelTimes[i4+2];
        auto tt  = mTravelTimes[i4+3];
        auto temtv = te - tv;
        // 1d operator (refracted time)
        auto t12min_1d = sycl::fmin(te + mDx*s2, tv + mDz*s1);
        // 2d operator
        auto t1_2d = mHuge;
        if (-temtv < mDx*s3 && temtv < mDz*s3)
        {
            auto sref2 = s3*s3;
            auto ta = tev + temtv;
            auto tb = tev - temtv;
            auto tab = ta - tb;
            auto tab2 = tab*tab;
            auto fourSref2 = four*sref2;
            t1_2d = (ta*mDx2i + tb*mDz2i)
                  + sycl::sqrt(fourSref2*mDx2i_p_Dz2i - mDx2i_Dz2i*tab2);
            t1_2d = t1_2d*mDx2i_p_Dz2i_inv;
        }
        mTravelTimesUpdate[i] = fmin(t12min_1d, t1_2d);
    }
private:
    const T *mTravelTimes;
    const T *mLocalSlowness;
    T *mTravelTimesUpdate;
    T mDx;
    T mDz;
    T mDx2i;
    T mDz2i;
    T mDx2i_Dz2i;
    T mDx2i_p_Dz2i;
    T mDx2i_p_Dz2i_inv;
    const T mHuge = std::numeric_limits<T>::max()/10;
    const T one = 1;
    const T four = 4;
};
*/

}
#endif
