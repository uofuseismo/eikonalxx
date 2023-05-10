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
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/io/vtkRectilinearGrid2d.hpp"
#include "private/grid.hpp"

using namespace EikonalXX;

namespace
{

template<typename T>
T getMin(const int n, const T x[])
{
    T vmin = 0;
//#ifdef ONE_API
//    sycl::queue queue{sycl::cpu_selector{}};
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
std::vector<T> interpolate2d(const size_t nGridX, const size_t nGridZ,
                             const T *velIn,
                             const EikonalXX::Ordering2D ordering)
{
    std::vector<T> slow;
    sycl::queue q{sycl::cpu_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> ();
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    // Figure out sizes and allocate space 
    auto nGrid = nGridX*nGridZ;
    auto nCellX = nGridX - 1;
    auto nCellZ = nGridZ - 1; 
    auto nCell = nCellX*nCellZ;
    auto velDevice = sycl::malloc_device<T> (nGrid, q);
    auto slowDevice =  sycl::malloc_device<T> (nCell, q);
    // Copy input velocities to device
    auto eCopyVelocity = q.submit([&](sycl::handler &h)
    {
        h.memcpy(velDevice, velIn, nGrid*sizeof(T)); 
    });
    // Interpolate to cell center.  The cell center is equidistant from all
    // nodes of the cell so a simple average is all that is required.
    sycl::range global{nCellX, nCellZ};
    auto nTileX = std::min(workGroupSize, nCellX);
    auto nTileZ = std::min(workGroupSize, nCellZ);
    sycl::range local{nTileX, nTileZ};
    const T four = 4;
    // Do interpolation
    auto eInterpolate = q.submit([&](sycl::handler &h)
    {
        h.depends_on(eCopyVelocity);
        if (ordering == EikonalXX::Ordering2D::NATURAL)
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<2> it)
            {
                size_t ix = it.get_global_id(0);
                size_t iz = it.get_global_id(1);
                auto idst = gridToIndex(nCellX, ix, iz); //iz*nCellX + ix;
                auto idx1 = iz*(nCellX + 1) + ix;
                auto idx2 = idx1 + 1;
                auto idx3 = idx1 + (nCellX + 1);
                auto idx4 = idx3 + 1;
                auto sum1 = velDevice[idx1] + velDevice[idx2];
                auto sum2 = velDevice[idx3] + velDevice[idx4];
                // 1/averageVelocity = 1/( (sum1 + sum2)/4 ) = 4/(sum1 + sum2)
                slowDevice[idst] = four/(sum1 + sum2);
            });
        }
        else
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<2> it)
            {
                size_t ix = it.get_global_id(0);
                size_t iz = it.get_global_id(1);
                auto idst = gridToIndex(nCellX, ix, iz); //iz*nCellX + ix;
                auto idx1 = ix*(nCellZ + 1) + iz;
                auto idx2 = idx1 + 1;
                auto idx3 = idx1 + (nCellZ + 1);
                auto idx4 = idx3 + 1;
                auto sum1 = velDevice[idx1] + velDevice[idx2];
                auto sum2 = velDevice[idx3] + velDevice[idx4];
                slowDevice[idst] = four/(sum1 + sum2);
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
    free(velDevice, q);
    free(slowDevice, q);
    return slow;
}

}

template<class T>
class Model2D<T>::Model2DImpl
{
public:
    /// Slowness field in each cell in seconds/meter
    Geometry2D mGeometry;
    std::vector<T> mSlowness;
    bool mHaveSlowness = false;
    bool mInitialized = false;
};

/// C'tor
template<class T>
Model2D<T>::Model2D() :
    pImpl(std::make_unique<Model2DImpl> ())
{
}

/// Copy c'tor
template<class T>
Model2D<T>::Model2D(const Model2D &model)
{
    *this = model;
}

/// Move c'tor
template<class T>
Model2D<T>::Model2D(Model2D &&model) noexcept
{
    *this = std::move(model);
}

/// Copy assignment
template<class T>
Model2D<T>& Model2D<T>::operator=(const Model2D<T> &model)
{
    if (&model == this){return *this;}
    pImpl = std::make_unique<Model2DImpl> (*model.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Model2D<T>& Model2D<T>::operator=(Model2D<T> &&model) noexcept
{
    if (&model == this){return *this;}
    pImpl = std::move(model.pImpl);
    return *this;
}

/// Destructor
template<class T>
Model2D<T>::~Model2D() = default;

/// Clear the class
template<class T>
void Model2D<T>::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mSlowness.clear();
    pImpl->mHaveSlowness = false;
    pImpl->mInitialized = false;
}

/// Initialize
template<class T>
void Model2D<T>::initialize(const Geometry2D &geometry)
{
    clear();
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Number of grid points in x not set");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Number of grid points in z not set");
    }
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x not set");
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
bool Model2D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Return geometry
template<class T>
Geometry2D Model2D<T>::getGeometry() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry;
}

/// Get the number of grid points
template<class T>
int Model2D<T>::getNumberOfGridPoints() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry.getNumberOfGridPoints();
}

/// Get the number of cells
template<class T>
int Model2D<T>::getNumberOfCells() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry.getNumberOfCells();
}

/// Set the nodal model -> this requires linear interpolation
template<class T>
//template<typename U>
void Model2D<T>::setNodalVelocities(const int nGridIn,
                                    const T v[],
                                    const EikonalXX::Ordering2D ordering)
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
    auto nz = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInZ());
    // Ensure all velocities are positive
    auto vmin = getMin(nx*nz, v);
    if (vmin <= 0)
    {
         throw std::invalid_argument("Non-positive velocities detected");
    }
    pImpl->mSlowness = interpolate2d(nx, nz, v, ordering);
#ifndef NDEBUG
    auto smin = getMin(pImpl->mSlowness.size(), pImpl->mSlowness.data());
    assert(smin > 0);
#endif
    pImpl->mHaveSlowness = true;
}

template<class T>
template<typename U>
void Model2D<T>::setCellularVelocities(const int nCellIn,
                                       const U v[],
                                       const EikonalXX::Ordering2D ordering)
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
    auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
    auto vmin = getMin(nCell, v);
    if (vmin <= 0)
    {
        throw std::invalid_argument("Non-positive velocities detected");
    }

    pImpl->mSlowness.resize(nCell, 0);
    const T one = 1;
    if (ordering == EikonalXX::Ordering2D::NATURAL)
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
        #pragma omp simd collapse(2)
        for (int iz=0; iz<nCellZ; ++iz)
        {
            for (int ix=0; ix<nCellX; ++ix)
            {
                auto isrc = gridToIndex(nCellZ, iz, ix); //ix*nCellZ + iz; 
                auto idst = gridToIndex(nCellX, ix, iz); //iz*nCellX + ix; 
                slow[idst] = one/v[isrc]; 
            }
        }
    }
    pImpl->mHaveSlowness = true;
}

/// Have slowness?
template<class T>
bool Model2D<T>::haveVelocities() const noexcept
{
    return pImpl->mHaveSlowness;
}

/// Get slowness model 
template<class T>
std::vector<T> Model2D<T>::getSlowness() const
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
const T* Model2D<T>::getSlownessPointer() const
{
    if (!haveVelocities())
    {
        throw std::runtime_error("Velocity model not yet set");
    }
    return pImpl->mSlowness.data();
}

/// Get the velocity model
template<class T>
std::vector<T> Model2D<T>::getVelocities() const
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
void Model2D<T>::writeVTK(const std::string &fileName,
                          const std::string &title) const
{
    auto v = getVelocities(); // Throws
    IO::VTKRectilinearGrid2D vtkWriter;
    constexpr bool writeBinary = true;
    vtkWriter.open(fileName, pImpl->mGeometry, title, writeBinary); 
    vtkWriter.writeCellularDataset(title, v.data(),
                                   EikonalXX::Ordering2D::NATURAL);
    vtkWriter.close();
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class EikonalXX::Model2D<double>;
template class EikonalXX::Model2D<float>;

//template void Eikonal::Model2D<double>::setNodalVelocities(
//    int nx, int nz, const double v[], Ordering2D ordering);
//template void Eikonal::Model2D<double>::setNodalVelocities(
//    int nx, int nz, const float v[], Ordering2D ordering);
template void EikonalXX::Model2D<double>::setCellularVelocities(
    int nCell, const double v[], Ordering2D ordering);
template void EikonalXX::Model2D<double>::setCellularVelocities(
    int nCell, const float v[], Ordering2D ordering);

//template void Eikonal::Model2D<float>::setNodalVelocities(
//    int nx, int nz, const double v[], Ordering2D ordering);
//template void Eikonal::Model2D<float>::setNodalVelocities(
//    int nx, int nz, const float v[], Ordering2D ordering);
template void EikonalXX::Model2D<float>::setCellularVelocities(
    int nCell, const double v[], Ordering2D ordering);
template void EikonalXX::Model2D<float>::setCellularVelocities(
    int nCell, const float v[], Ordering2D ordering);
