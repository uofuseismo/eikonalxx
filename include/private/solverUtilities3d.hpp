#ifndef PRIVATE_SOLVERUTILITIES3D_HPP
#define PRIVATE_SOLVERUTILITIES3D_HPP
#include <CL/sycl.hpp>
#include <cmath>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#include <tbb/tbb.h>
#include "eikonalxx/enums.hpp"
#include "private/grid.hpp"
#include "private/pad.hpp"

#define BOUNDARY_NODE -1
#define SOURCE_NODE 0
#define UPDATE_NODE 1
namespace
{

using namespace EikonalXX;

/// @brief Container for the slownesses in a level for the 2D solver.
template<class T>
struct SweepSlowness3D
{
    /// C'tor
    SweepSlowness3D() = default;
    /// Copy c'tor
    SweepSlowness3D(const SweepSlowness3D &slowness)
    {
        *this = slowness;
    }
    /// Move c'tor
    SweepSlowness3D(SweepSlowness3D &&slowness) noexcept
    {
        *this = std::move(slowness);
    }
    /// Copy assingment
    SweepSlowness3D& operator=(const SweepSlowness3D &slowness)
    {
        if (&slowness == this){return *this;}
        if (slowness.nNodes > 0)
        {
            allocate(nNodes);
            std::copy(slowness.s0, s0, nNodes);
            std::copy(slowness.s1, s1, nNodes);
            std::copy(slowness.s3, s3, nNodes);
            std::copy(slowness.s4, s4, nNodes);
            std::copy(slowness.s5, s5, nNodes);
            std::copy(slowness.s6, s6, nNodes);
            std::copy(slowness.s7, s7, nNodes);
        }
        else
        {
            nNodes = 0;
            s0 = nullptr;
            s1 = nullptr;
            s3 = nullptr;
            s4 = nullptr;
            s5 = nullptr;
            s6 = nullptr;
            s7 = nullptr;
        }
        return *this;
    }
    /// Move assignment
    SweepSlowness3D& operator=(SweepSlowness3D &&slowness) noexcept
    {
        if (&slowness == this){return *this;}
        nNodes = slowness.nNodes;
        s0 = std::move(slowness.s0);
        s1 = std::move(slowness.s1);
        s3 = std::move(slowness.s3);
        s4 = std::move(slowness.s4);
        s5 = std::move(slowness.s5);
        s6 = std::move(slowness.s6);
        s7 = std::move(slowness.s7);
        return *this;
    }
    /// Destructor
    ~SweepSlowness3D()
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
        auto nBytes = static_cast<size_t> (nPad)*sizeof(T);
        s0 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s1 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s3 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s4 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s5 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s6 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        s7 = static_cast<T *> (std::aligned_alloc(64, nBytes));
        std::fill(s0, s0 + nPad, 0);
        std::fill(s1, s1 + nPad, 0);
        std::fill(s3, s3 + nPad, 0);
        std::fill(s4, s4 + nPad, 0);
        std::fill(s5, s5 + nPad, 0);
        std::fill(s6, s6 + nPad, 0);
        std::fill(s7, s7 + nPad, 0);
    }    
    /// This sets all the slownesses in the level to 0
    void zero() noexcept
    {
        if (nNodes > 0) 
        {
            std::fill(s0, s0 + nNodes, 0);
            std::fill(s1, s1 + nNodes, 0);
            std::fill(s3, s3 + nNodes, 0);
            std::fill(s4, s4 + nNodes, 0);
            std::fill(s5, s5 + nNodes, 0);
            std::fill(s6, s6 + nNodes, 0);
            std::fill(s7, s7 + nNodes, 0);
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
            sMin = std::min(sMin, *std::min_element(s4, s4 + nNodes));
            sMin = std::min(sMin, *std::min_element(s5, s5 + nNodes));
            sMin = std::min(sMin, *std::min_element(s6, s6 + nNodes));
            sMin = std::min(sMin, *std::min_element(s7, s7 + nNodes));
        }
        return sMin;
    }    
    /// Releases memory.
    void clear() noexcept
    {
        if (s0){free(s0);}
        if (s1){free(s1);}
        if (s3){free(s3);}
        if (s4){free(s4);}
        if (s5){free(s5);}
        if (s6){free(s6);}
        if (s7){free(s7);}
        s0 = nullptr;
        s1 = nullptr;
        s3 = nullptr;
        s4 = nullptr;
        s5 = nullptr;
        s6 = nullptr;
        s7 = nullptr;
        nNodes = 0; 
    }    
    /// The slowness (s/m) in the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s0 = nullptr;
    /// The slowness (s/m) in the cell to the right or left of the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s1 = nullptr;
    /// The slowness (s/m) in the cell in front of or in back of the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s3 = nullptr;
    /// The slowness (s/m) in the cell above or below the home cell.
    /// This is an array whose dimension is at least [nNodes].
    T *__attribute__((aligned(64))) s4 = nullptr;
    /// The slowness (s/m) in the cell above or below cell s1.
    T *__attribute__((aligned(64))) s5 = nullptr;
    /// The slowness (s/m) in the cell above or below cell s3.
    T *__attribute__((aligned(64))) s6 = nullptr;
    T *__attribute__((aligned(64))) s7 = nullptr;
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
void getSweepFiniteDifferenceSigns(//const EikonalXX::SweepNumber3D sweep,
                                   int *ixShift, int *iyShift, int *izShift,
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
                 int *level)//, int *indx)
{
    int jx, jy, jz;
    permuteGrid<E>(ix, iy, iz, nx, ny, nz, &jx, &jy, &jz);
    //*indx = jz;
    *level = jx + jy + jz;
}


}
#endif
