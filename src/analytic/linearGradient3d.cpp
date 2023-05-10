#include <string>
#include <CL/sycl.hpp>
#include "eikonalxx/analytic/linearGradient3d.hpp"
#include "eikonalxx/source3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/io/vtkRectilinearGrid3d.hpp"
#include "private/grid.hpp"

namespace
{

template<typename T>
T acoshArgument(const T delX, const T delY, const T delZ,
                const T vel0, const T slow0, const T absG2,
                const T vGradInX, const T vGradInY, const T vGradInZ)
{
    T dist2 = delX*delX + delY*delY + delZ*delZ;
    T vel = vel0
          + vGradInX*delX
          + vGradInY*delY
          + vGradInZ*delZ;
    T arg = 1 + (slow0*absG2*dist2)/(2*vel);
    return arg;
}

/// Solves the eikonal equation in a 3D linear-gradient
template<class T>
void solve3d(const size_t nGridX, const size_t nGridY, const size_t nGridZ,
             const double xSrcOffset,
             const double ySrcOffset,
             const double zSrcOffset,
             const double dx, const double dy, const double dz,
             const double velTop, const double velBottom,
             std::vector<T> *travelTimes)
{
    sycl::queue q{sycl::cpu_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> (); 
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    constexpr T vGradInX = 0;
    constexpr T vGradInY = 0;
    T vGradInZ = static_cast<T> ((velBottom - velTop)/((nGridZ - 1)*dz));
    T absG2 = vGradInX*vGradInX + vGradInY*vGradInY + vGradInZ*vGradInZ; 
    T absG = std::sqrt(absG2);
    T vel0 = static_cast<T> (velTop + vGradInZ*zSrcOffset);
    auto slow0 = static_cast<T> (1/vel0);
    // Figure out sizes and allocate space 
    auto nGrid = nGridX*nGridY*nGridZ;
    auto travelTimesDevice = sycl::malloc_device<T> (nGrid, q);
    // Determine ranges (with tiling)
    sycl::range global{nGridX, nGridY, nGridZ};
    auto nTileX = std::min(workGroupSize, nGridX);
    auto nTileY = std::min(workGroupSize, nGridY);
    auto nTileZ = std::min(workGroupSize, nGridZ);
    sycl::range local{nTileX, nTileY, nTileZ};
    // Simplify some geometric terms
    auto xShiftedSource = static_cast<T> (xSrcOffset);
    auto yShiftedSource = static_cast<T> (ySrcOffset);
    auto zShiftedSource = static_cast<T> (zSrcOffset);
    auto hx = static_cast<T> (dx);
    auto hy = static_cast<T> (dy);
    auto hz = static_cast<T> (dz);
    // Compute L2 distance from the source to each point in grid
    // and scale by slowness.
    auto eTravelTimes = q.submit([&](sycl::handler &h) 
    {
        if (std::abs(velTop - velBottom) > 1.e-14)
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<3> it) 
            {
                size_t ix = it.get_global_id(0);
                size_t iy = it.get_global_id(1);
                size_t iz = it.get_global_id(2);
                auto idst = ::gridToIndex(nGridX, nGridY, ix, iy, iz);
                T delX = ix*hx - xShiftedSource;
                T delY = iy*hy - yShiftedSource;
                T delZ = iz*hz - zShiftedSource;
                T arg = ::acoshArgument(delX, delY, delZ,
                                        vel0, slow0, absG2,
                                        vGradInX, vGradInY, vGradInZ);
                travelTimesDevice[idst] = sycl::acosh(arg)/absG;
            });
        }
        else
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<3> it) 
            {   
                size_t ix = it.get_global_id(0);
                size_t iy = it.get_global_id(1);
                size_t iz = it.get_global_id(2);
                auto idst = ::gridToIndex(nGridX, nGridY, ix, iy, iz);
                T delX = xShiftedSource - ix*hx;
                T delY = yShiftedSource - iy*hy;
                T delZ = zShiftedSource - iz*hz;
                travelTimesDevice[idst]
                      = sycl::sqrt(delX*delX + delY*delY + delZ*delZ)*slow0;
            });
        }
    });
    // Return slowness field to host
    if (travelTimes->size() != nGrid)
    {
        travelTimes->resize(nGrid, 0);
    }
    auto travelTimesPtr = travelTimes->data();
    q.submit([&](sycl::handler &h) 
    {
        h.depends_on(eTravelTimes);
        h.memcpy(travelTimesPtr, travelTimesDevice, nGrid*sizeof(T));
    }); 
    q.wait();
    // Release memory
    sycl::free(travelTimesDevice, q);
}
}

using namespace EikonalXX::Analytic;

template<class T>
class LinearGradient3D<T>::LinearGradient3DImpl
{
public:
    EikonalXX::Geometry3D mGeometry;
    EikonalXX::Source3D mSource;
    std::vector<T> mTravelTimeField; 
    std::pair<double, double> mVelocity{0, 0};
    bool mHaveTravelTimeField{false};
    bool mHaveSource{false};
    bool mInitialized{false};
};


/// C'tor
template<class T>
LinearGradient3D<T>::LinearGradient3D() :
    pImpl(std::make_unique<LinearGradient3DImpl> ())
{
}

/// Copy c'tor
template<class T>
LinearGradient3D<T>::LinearGradient3D(const LinearGradient3D &solver)
{
    *this = solver;
}

/// Move c'tor
template<class T>
LinearGradient3D<T>::LinearGradient3D(LinearGradient3D &&solver) noexcept
{
    *this = std::move(solver);
}

/// Reset the class
template<class T>
void LinearGradient3D<T>::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mSource.clear();
    pImpl->mTravelTimeField.clear();
    pImpl->mVelocity = std::pair<double, double> {0, 0};
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveSource = false;
    pImpl->mInitialized = false;
}

/// Destructor
template<class T>
LinearGradient3D<T>::~LinearGradient3D() = default;

/// Copy assignment
template<class T>
LinearGradient3D<T>& LinearGradient3D<T>::operator=(
    const LinearGradient3D &solver)
{
    if (&solver == this){return *this;}
    pImpl = std::make_unique<LinearGradient3DImpl> (*solver.pImpl);
    return *this;
}

/// Move assignment
template<class T>
LinearGradient3D<T>& LinearGradient3D<T>::operator=(
    LinearGradient3D &&solver) noexcept
{
    if (&solver == this){return *this;}
    pImpl = std::move(solver.pImpl);
    return *this;
}

/// Initialize the class
template<class T>
void LinearGradient3D<T>::initialize(const Geometry3D &geometry)
{
    clear();
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Grid points in x not set on geometry");
    }
    if (!geometry.haveNumberOfGridPointsInY())
    {
        throw std::invalid_argument("Grid points in y not set on geometry");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Grid points in z not set on geometry");
    }
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x not set on geometry");
    }
    if (!geometry.haveGridSpacingInY())
    {
        throw std::invalid_argument("Grid spacing in y not set on geometry");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("Grid spacing in z not set on geometry");
    }
    pImpl->mGeometry = geometry;
    auto nGrid = static_cast<size_t> (geometry.getNumberOfGridPoints());
    pImpl->mTravelTimeField.resize(nGrid, 0);
    pImpl->mInitialized = true;
}

/// Initialized?
template<class T>
bool LinearGradient3D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Set source
template<class T>
void LinearGradient3D<T>::setSource(const EikonalXX::Source3D &source)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!source.haveLocationInX())
    {
        throw std::invalid_argument("Source location in x not set");
    }
    if (!source.haveLocationInY())
    {
        throw std::invalid_argument("Source location in y not set");
    }
    if (!source.haveLocationInZ())
    {
        throw std::invalid_argument("Source location in z not set");
    }
    pImpl->mSource = source;
    pImpl->mHaveSource = true;
}

/// Set the source from an (x,y, z) pair
template<class T>
void LinearGradient3D<T>::setSource(
    const std::tuple<double, double, double> &sourceLocation)
{
    pImpl->mHaveSource = false; 
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dy = pImpl->mGeometry.getGridSpacingInY();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    auto ny = pImpl->mGeometry.getNumberOfGridPointsInY();
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
    auto xmin = pImpl->mGeometry.getOriginInX();
    auto ymin = pImpl->mGeometry.getOriginInY();
    auto zmin = pImpl->mGeometry.getOriginInZ();
    auto xmax = xmin + (nx - 1)*dx;
    auto ymax = ymin + (ny - 1)*dy;
    auto zmax = zmin + (nz - 1)*dz; 
    if (std::get<0> (sourceLocation) < xmin ||
        std::get<0> (sourceLocation) > xmax)
    {
        throw std::invalid_argument("x source position = "
                                  + std::to_string(std::get<0> (sourceLocation))
                                  + " must be in range [" + std::to_string(xmin)
                                  + "," + std::to_string(xmax) + "]");
    }
    if (std::get<1> (sourceLocation) < ymin ||
        std::get<1> (sourceLocation) > ymax)
    {
        throw std::invalid_argument("y source position = "
                                  + std::to_string(std::get<1> (sourceLocation))
                                  + " must be in range [" + std::to_string(ymin)
                                  + "," + std::to_string(ymax) + "]");
    }
    if (std::get<2> (sourceLocation) < zmin ||
        std::get<2> (sourceLocation) > zmax)
    {
        throw std::invalid_argument("z source position = "
                                  + std::to_string(std::get<2> (sourceLocation))
                                  + " must be in range [" + std::to_string(zmin)
                                  + "," + std::to_string(zmax) + "]");
    }
    // Create a source
    EikonalXX::Source3D source;
    source.setGeometry(pImpl->mGeometry);
    source.setLocationInX(std::get<0> (sourceLocation));
    source.setLocationInY(std::get<1> (sourceLocation));
    source.setLocationInZ(std::get<2> (sourceLocation));
    setSource(source);
}

/// Have source?
template<class T>
bool LinearGradient3D<T>::haveSource() const noexcept
{
    return pImpl->mHaveSource;
}

/// Get source
template<class T>
EikonalXX::Source3D LinearGradient3D<T>::getSource() const
{
    if (!haveSource())
    {
        throw std::runtime_error("Source location not set");
    }
    return pImpl->mSource;
}

/// Set velocity model
template<class T>
void LinearGradient3D<T>::setVelocityModel(
    const std::pair<double, double> &velocity)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (velocity.first <= 0)
    {
        throw std::invalid_argument("velocity.first = "
                                  + std::to_string(velocity.first)
                                  + " must be positive");
    }
    if (velocity.second <= 0)
    {
        throw std::invalid_argument("velocity.second = "
                                  + std::to_string(velocity.second)
                                  + " must be positive");
    }
    pImpl->mVelocity = velocity;
}

/// Have velocity?
template<class T>
bool LinearGradient3D<T>::haveVelocityModel() const noexcept
{
    return (pImpl->mVelocity.first > 0 && pImpl->mVelocity.second > 0);
}

/// Solve the eikonal equation
template<class T>
void LinearGradient3D<T>::solve()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!haveSource()){throw std::runtime_error("Source not yet set");}
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    pImpl->mHaveTravelTimeField = false;
    auto nx = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInX());
    auto ny = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInY());
    auto nz = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInZ());
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dy = pImpl->mGeometry.getGridSpacingInY();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto xSrcOffset = pImpl->mSource.getOffsetInX();
    auto ySrcOffset = pImpl->mSource.getOffsetInY();
    auto zSrcOffset = pImpl->mSource.getOffsetInZ();
    // Solve it
    ::solve3d(nx, ny, nz,
              xSrcOffset, ySrcOffset, zSrcOffset,
              dx, dy, dz,
              pImpl->mVelocity.first, pImpl->mVelocity.second,
              &pImpl->mTravelTimeField);
    pImpl->mHaveTravelTimeField = true;
}

/// Have travel time field?
template<class T>
bool LinearGradient3D<T>::haveTravelTimeField() const noexcept
{
    return pImpl->mHaveTravelTimeField;
}

template<class T>
std::vector<T> LinearGradient3D<T>::getTravelTimeField() const
{
    if (!haveTravelTimeField())
    {
        throw std::runtime_error("Travel time field not yet computed");
    }
    return pImpl->mTravelTimeField;
}

template<class T>
const T* LinearGradient3D<T>::getTravelTimeFieldPointer() const
{
    if (!haveTravelTimeField())
    {
        throw std::runtime_error("Travel time field not yet computed");
    }
    return pImpl->mTravelTimeField.data();
}

/// Write the travel time field
template<class T>
void LinearGradient3D<T>::writeVTK(const std::string &fileName,
                                   const std::string &title) const
{
    auto tPtr = getTravelTimeFieldPointer(); // Throws
    IO::VTKRectilinearGrid3D vtkWriter;
    constexpr bool writeBinary = true;
    vtkWriter.open(fileName, pImpl->mGeometry, title, writeBinary); 
    vtkWriter.writeNodalDataset(title, tPtr, EikonalXX::Ordering3D::Natural);
    vtkWriter.close();
}

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///

template class EikonalXX::Analytic::LinearGradient3D<double>;
template class EikonalXX::Analytic::LinearGradient3D<float>;
