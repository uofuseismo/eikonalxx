#include <CL/sycl.hpp>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/geometry2d.hpp"

namespace
{
/// Solves the eikonal equation in a 2D homogeneous medium.
template<class T>
std::vector<T> solve2d(const size_t nGridX, const size_t nGridZ,
                       const double xSrc, const double zSrc,
                       const double x0, const double z0,
                       const double dx, const double dz,
                       const double vel)
{
    std::vector<T> ttimes;
    sycl::queue q{sycl::cpu_selector{},
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> (); 
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    // Figure out sizes and allocate space 
    auto nGrid = nGridX*nGridZ;
    auto ttimesDevice = sycl::malloc_device<T> (nGrid, q); 
    // Determine ranges (with tiling)
    sycl::range global{nGridZ, nGridX};
    auto nTileX = std::min(workGroupSize, nGridX);
    auto nTileZ = std::min(workGroupSize, nGridZ);
    sycl::range local{nTileX, nTileZ};
    // Simplify some geometric terms
    auto xShiftedSource = static_cast<T> (xSrc - x0);
    auto zShiftedSource = static_cast<T> (zSrc - z0);
    auto hx = static_cast<T> (dx);
    auto hz = static_cast<T> (dz);
    auto slowness = static_cast<T> (1./vel); // Do division here
    // Compute L2 distance from the source to each point in grid
    // and scale by slowness.
    auto eTravelTimes = q.submit([&](sycl::handler &h) 
    {
        h.parallel_for(sycl::nd_range{global, local},
                       [=](sycl::nd_item<2> it) 
        {
            size_t iz = it.get_global_id(0);
            size_t ix = it.get_global_id(1);
            auto idst = iz*nGridX + ix;
            T delX = xShiftedSource - ix*hx;
            T delZ = zShiftedSource - iz*hz;
            ttimesDevice[idst] = sycl::hypot(delX, delZ)*slowness;
        });
    });
    // Return slowness field to host
    ttimes.resize(nGrid, 0);
    auto ttimesPtr = ttimes.data();
    q.submit([&](sycl::handler &h) 
    {
        h.depends_on(eTravelTimes);
        h.memcpy(ttimesPtr, ttimesDevice, nGrid*sizeof(T));
    }); 
    q.wait();
    // Release memory
    free(ttimesDevice, q); 
    return ttimes;
}
}

using namespace EikonalXX::Analytic;

template<class T>
class Homogeneous2D<T>::Homogeneous2DImpl
{
public:
    EikonalXX::Geometry2D mGeometry;
    EikonalXX::Source2D mSource;
    std::vector<T> mTravelTimeField; 
    bool mInitialized = false; 
};


/// C'tor
template<class T>
Homogeneous2D<T>::Homogeneous2D() :
    pImpl(std::make_unique<Homogeneous2DImpl> ())
{
}

/// Copy c'tor
template<class T>
Homogeneous2D<T>::Homogeneous2D(const Homogeneous2D &solver)
{
    *this = solver;
}

/// Move c'tor
template<class T>
Homogeneous2D<T>::Homogeneous2D(Homogeneous2D &&solver) noexcept
{
    *this = std::move(solver);
}

/// Destructor
template<class T>
Homogeneous2D<T>::~Homogeneous2D() = default;

/// Copy assignment
template<class T>
Homogeneous2D<T>& Homogeneous2D<T>::operator=(const Homogeneous2D &solver)
{
    if (&solver == this){return *this;}
    pImpl = std::make_unique<Homogeneous2DImpl> (*solver.pImpl);
    return *this;
}

/// Move assignment
template<class T>
Homogeneous2D<T>& Homogeneous2D<T>::operator=(Homogeneous2D &&solver) noexcept
{
    if (&solver == this){return *this;}
    pImpl = std::move(solver.pImpl);
    return *this;
}

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///

template class EikonalXX::Analytic::Homogeneous2D<double>;
template class EikonalXX::Analytic::Homogeneous2D<float>;
