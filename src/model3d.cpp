#include <iostream>
#include <CL/sycl.hpp>
#include <algorithm>
#include <cassert>
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
#include "eikonalxx/model3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/io/vtkRectilinearGrid3d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX;

namespace
{

template<typename T>
T getMin(const int n, const T x[])
{
    T vmin = 0;
//#ifdef ONE_API
//    sycl::queue queue{sycl::cpu_selector_v};
//    auto policy = dpstd::execution::make_device_policy<class Min>(queue);
//    vmin = *std::min_element(policy, x, x + n);
//    queue.wait();
//#else
#ifdef USE_PSTL
    vmin = *std::min_element(std::execution::unseq, x, x + n);
#else
    vmin = *std::min_element(x, x + n);
#endif
    return vmin;
}

/// Interpolates a nodal velocity model onto a cell-based slowness model
template<class T>
std::vector<T> interpolate3d(const size_t nGridX,
                             const size_t nGridY,
                             const size_t nGridZ,
                             const T *velIn,
                             const EikonalXX::Ordering3D ordering)
{
    std::vector<T> slow;
    sycl::queue q{sycl::cpu_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> ();
    workGroupSize = static_cast<size_t> (std::cbrt(workGroupSize));
    // Figure out sizes and allocate space 
    auto nGrid = nGridX*nGridY*nGridZ;
    auto nCellX = nGridX - 1;
    auto nCellY = nGridY - 1;
    auto nCellZ = nGridZ - 1;
    auto nCell = nCellX*nCellY*nCellZ;
    auto velDevice = sycl::malloc_device<T> (nGrid, q);
    auto slowDevice =  sycl::malloc_device<T> (nCell, q);
    // Copy input velocities to device
    auto eCopyVelocity = q.submit([&](sycl::handler &h)
    {
        h.memcpy(velDevice, velIn, nGrid*sizeof(T)); 
    });
    // Interpolate to cell center.  The cell center is equidistant from all
    // nodes of the cell so a simple average is all that is required.
    sycl::range global{nCellX, nCellY, nCellZ};
    auto nTileX = std::min(workGroupSize, nCellX);
    auto nTileY = std::min(workGroupSize, nCellY);
    auto nTileZ = std::min(workGroupSize, nCellZ);
    sycl::range local{nTileX, nTileY, nTileZ};
    const T eight = 8;
    // Do interpolation
    auto eInterpolate = q.submit([&](sycl::handler &h)
    {
        h.depends_on(eCopyVelocity);
        if (ordering == EikonalXX::Ordering3D::Natural)
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<3> it)
            {
                size_t ix = it.get_global_id(0);
                size_t iy = it.get_global_id(1);
                size_t iz = it.get_global_id(2);
                auto idst = gridToIndex(nCellX, nCellY, ix, iy, iz); //iz*(nCellX*nCellY) + iy*nCellX + ix;
                auto idx1 = iz*(nCellX + 1)*(nCellY + 1)
                          + iy*(nCellX + 1)
                          + ix;
                auto idx2 = idx1 + 1;
                auto idx3 = idx1 + (nCellX + 1);
                auto idx4 = idx3 + 1;
                auto idx5 = idx1 + (nCellX + 1)*(nCellY + 1);
                auto idx6 = idx2 + (nCellX + 1)*(nCellY + 1);
                auto idx7 = idx3 + (nCellX + 1)*(nCellY + 1);
                auto idx8 = idx4 + (nCellX + 1)*(nCellY + 1);
                auto sum1 = velDevice[idx1] + velDevice[idx2];
                auto sum2 = velDevice[idx3] + velDevice[idx4];
                auto sum3 = velDevice[idx5] + velDevice[idx6];
                auto sum4 = velDevice[idx7] + velDevice[idx8];
                // 1/averageVelocity = 1/( (sum1 + sum2 + sum3 + sum4)/8 )
                //                   = 8/(sum1 + sum2 + sum3 + sum4)
                slowDevice[idst] = eight/(sum1 + sum2 + sum3 + sum4);
            });
        }
        else
        {
            h.depends_on(eCopyVelocity);
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<3> it)
            {
                size_t ix = it.get_global_id(0);
                size_t iy = it.get_global_id(1);
                size_t iz = it.get_global_id(2);
                auto idst = gridToIndex(nCellX, nCellY, ix, iy, iz); //iz*(nCellX*nCellY) + iy*nCellX + ix;
                auto idx1 = ix*(nCellZ + 1)*(nCellY + 1)
                          + iy*(nCellZ + 1)
                          + iz;
                auto idx2 = idx1 + 1;
                auto idx3 = idx1 + (nCellZ + 1);
                auto idx4 = idx3 + 1;
                auto idx5 = idx1 + (nCellZ + 1)*(nCellY + 1);
                auto idx6 = idx2 + (nCellZ + 1)*(nCellY + 1);
                auto idx7 = idx3 + (nCellZ + 1)*(nCellY + 1);
                auto idx8 = idx4 + (nCellZ + 1)*(nCellY + 1);
                auto sum1 = velDevice[idx1] + velDevice[idx2];
                auto sum2 = velDevice[idx3] + velDevice[idx4];
                auto sum3 = velDevice[idx5] + velDevice[idx6];
                auto sum4 = velDevice[idx7] + velDevice[idx8];
                slowDevice[idst] = eight/(sum1 + sum2 + sum3 + sum4);
            });
        }
    });
    // Return slowness field to host
    slow.resize(nCell, -1);
    auto slowPtr = slow.data();
    q.submit([&](sycl::handler &h)
    {
        h.depends_on(eInterpolate);
        h.memcpy(slowPtr, slowDevice, nCell*sizeof(T));
    });
    q.wait();
    // Release memory
    sycl::free(velDevice, q);
    sycl::free(slowDevice, q);
    return slow;
}

}

template<class T>
class Model3D<T>::Model3DImpl
{
public:
    /// Slowness field in each cell in seconds/meter
    Geometry3D mGeometry;
    std::vector<T> mSlowness;
    bool mHaveSlowness = false;
    bool mInitialized = false;
};

/// C'tor
template<class T>
Model3D<T>::Model3D() :
    pImpl(std::make_unique<Model3DImpl> ())
{
}

/// Copy c'tor
template<class T>
Model3D<T>::Model3D(const Model3D &model)
{
    *this = model;
}

/// Move c'tor
template<class T>
Model3D<T>::Model3D(Model3D &&model) noexcept
{
    *this = std::move(model);
}

/// Copy assignment
template<class T>
Model3D<T>& Model3D<T>::operator=(const Model3D<T> &model)
{
    if (&model == this){return *this;}
    pImpl = std::make_unique<Model3DImpl> (*model.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Model3D<T>& Model3D<T>::operator=(Model3D<T> &&model) noexcept
{
    if (&model == this){return *this;}
    pImpl = std::move(model.pImpl);
    return *this;
}

/// Destructor
template<class T>
Model3D<T>::~Model3D() = default;

/// Clear the class
template<class T>
void Model3D<T>::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mSlowness.clear();
    pImpl->mHaveSlowness = false;
    pImpl->mInitialized = false;
}

/// Initialize
template<class T>
void Model3D<T>::initialize(const Geometry3D &geometry)
{
    clear();
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Number of grid points in x not set");
    }
    if (!geometry.haveNumberOfGridPointsInY())
    {
        throw std::invalid_argument("Number of grid points in y not set");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Number of grid points in z not set");
    }
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x not set");
    }
    if (!geometry.haveGridSpacingInY())
    {   
        throw std::invalid_argument("Grid spacing in y not set");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("Grid spacing in z not set");
    }
    pImpl->mGeometry = geometry;
    pImpl->mInitialized = true; 
}

/// Is class initialized?
template<class T>
bool Model3D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Return geometry
template<class T>
Geometry3D Model3D<T>::getGeometry() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry;
}

/// Get the number of grid points
template<class T>
int Model3D<T>::getNumberOfGridPoints() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry.getNumberOfGridPoints();
}

/// Get the number of cells
template<class T>
int Model3D<T>::getNumberOfCells() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry.getNumberOfCells();
}

/// Set the nodal model -> this requires linear interpolation
template<class T>
//template<typename U>
void Model3D<T>::setNodalVelocities(const int nGridIn,
                                    const T v[],
                                    const EikonalXX::Ordering3D ordering)
{
    pImpl->mHaveSlowness = false;
    auto n = getNumberOfGridPoints(); // Throws
    if (nGridIn != n)
    {
        throw std::invalid_argument("Number of grid points = "
                                  + std::to_string(nGridIn)
                                  + " must equal " + std::to_string(n));
    }
    if (v == nullptr){throw std::invalid_argument("velocities is NULL");}
    auto nx = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInX());
    auto ny = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInY());
    auto nz = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInZ());
    // Ensure all velocities are positive
    auto vmin = getMin(nx*ny*nz, v);
    if (vmin <= 0)
    {
         throw std::invalid_argument("Non-positive velocities detected");
    }
    pImpl->mSlowness = interpolate3d(nx, ny, nz, v, ordering);
#ifndef NDEBUG
    auto smin = getMin(pImpl->mSlowness.size(), pImpl->mSlowness.data());
    assert(smin > 0);
#endif
    pImpl->mHaveSlowness = true;
}

template<class T>
template<typename U>
void Model3D<T>::setCellularVelocities(const int nCellIn,
                                       const U v[],
                                       const EikonalXX::Ordering3D ordering)
{
    pImpl->mHaveSlowness = false;
    auto nCell = getNumberOfCells();
    if (nCellIn != nCell)
    {
        throw std::invalid_argument("Number of cells = "
                                  + std::to_string(nCellIn)
                                  + " must equal " + std::to_string(nCell));
    }
    if (v == nullptr){throw std::invalid_argument("velocities is NULL");}
    // Check the velocities
    auto nCellX = pImpl->mGeometry.getNumberOfCellsInX();
    auto nCellY = pImpl->mGeometry.getNumberOfCellsInY();
    auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
    auto vmin = getMin(nCell, v);
    if (vmin <= 0)
    {
        throw std::invalid_argument("Non-positive velocities detected");
    }

    pImpl->mSlowness.resize(nCell, 0);
    const T one = 1;
    if (ordering == EikonalXX::Ordering3D::Natural)
    {
#ifdef USE_PSTL
        std::transform(std::execution::unseq,
                       v, v + nCell, pImpl->mSlowness.begin(),
                       std::bind1st(std::divides<U> (), one)); 
#else
        std::transform(v, v + nCell, pImpl->mSlowness.begin(),
                       std::bind1st(std::divides<U> (), one));
#endif
    }
    else
    {
        // Transpose and compute slowness 
        auto slow = pImpl->mSlowness.data();
        #pragma omp simd collapse(3)
        for (int iz = 0; iz < nCellZ; ++iz)
        {
            for (int iy = 0; iy < nCellY; ++iy)
            {
                for (int ix=0; ix < nCellX; ++ix)
                {
                    // ix*nCellZ*nCellY + iy*nCellZ + iz
                    auto isrc = gridToIndex(nCellZ, nCellY, iz, iy, ix);
                    // iz*nCellX*nCellY + iy*nCellX + ix
                    auto idst = gridToIndex(nCellX, nCellY, ix, iy, iz);
                    slow[idst] = one/v[isrc]; 
                }
            }
        }
    }
    pImpl->mHaveSlowness = true;
}

/// Have slowness?
template<class T>
bool Model3D<T>::haveVelocities() const noexcept
{
    return pImpl->mHaveSlowness;
}

/// Get slowness model 
template<class T>
std::vector<T> Model3D<T>::getSlowness() const
{
    auto slownessPtr = getSlownessPointer(); // Throws
    auto nCells = getNumberOfCells();
    std::vector<T> slowness(nCells);
#ifdef USE_PSTL
    std::copy(std::execution::unseq,
              slownessPtr, slownessPtr + nCells, slowness.begin());
#else
    std::copy(slownessPtr, slownessPtr + nCells, slowness.begin());
#endif 
    return slowness;
}

/// Get a pointer to the slowness model
template<class T>
const T* Model3D<T>::getSlownessPointer() const
{
    if (!haveVelocities())
    {
        throw std::runtime_error("Velocity model not yet set");
    }
    return pImpl->mSlowness.data();
}

/// Get the velocity model
template<class T>
std::vector<T> Model3D<T>::getVelocities() const
{
    if (!haveVelocities())
    {
        throw std::runtime_error("Velocity model not yet set");
    }
    const T one = 1;
    std::vector<T> velocities(pImpl->mSlowness.size());
#ifdef USE_PSTL
    std::transform(std::execution::unseq,
                   pImpl->mSlowness.begin(), pImpl->mSlowness.end(),
                   velocities.begin(),
                   std::bind1st(std::divides<T> (), one));
#else
    std::transform(pImpl->mSlowness.begin(), pImpl->mSlowness.end(),
                   velocities.begin(),
                   std::bind1st(std::divides<T> (), one));
#endif
    return velocities;
}

/// Write the velocity model
template<class T>
void Model3D<T>::writeVTK(const std::string &fileName,
                          const std::string &title) const
{
    auto v = getVelocities(); // Throws
    IO::VTKRectilinearGrid3D vtkWriter;
    constexpr bool writeBinary = true;
    vtkWriter.open(fileName, pImpl->mGeometry, title, writeBinary); 
    vtkWriter.writeCellularDataset(title, v.data(),
                                   EikonalXX::Ordering3D::Natural);
    vtkWriter.close();
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class EikonalXX::Model3D<double>;
template class EikonalXX::Model3D<float>;

template void EikonalXX::Model3D<double>::setCellularVelocities(
    int nCell, const double v[], Ordering3D ordering);
template void EikonalXX::Model3D<double>::setCellularVelocities(
    int nCell, const float v[], Ordering3D ordering);

template void EikonalXX::Model3D<float>::setCellularVelocities(
    int nCell, const double v[], Ordering3D ordering);
template void EikonalXX::Model3D<float>::setCellularVelocities(
    int nCell, const float v[], Ordering3D ordering);
