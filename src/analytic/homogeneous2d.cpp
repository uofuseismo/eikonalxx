#include <string>
#include <CL/sycl.hpp>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/io/vtkRectilinearGrid2d.hpp"
#include "private/grid.hpp"

namespace
{
/// Solves the eikonal equation in a 2D homogeneous medium.
template<typename T>
void solve2d(const size_t nGridX, const size_t nGridZ,
             const double xSrcOffset, const double zSrcOffset,
             const double dx, const double dz,
             const double vel,
             std::vector<T> *travelTimes)
{
    sycl::queue q{sycl::cpu_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> (); 
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    // Figure out sizes and allocate space 
    auto nGrid = nGridX*nGridZ;
    auto travelTimesDevice = sycl::malloc_device<T> (nGrid, q); 
    // Determine ranges (with tiling)
    sycl::range global{nGridX, nGridZ};
    auto nTileX = std::min(workGroupSize, nGridX);
    auto nTileZ = std::min(workGroupSize, nGridZ);
    sycl::range local{nTileX, nTileZ};
    // Simplify some geometric terms
    auto xShiftedSource = static_cast<T> (xSrcOffset);
    auto zShiftedSource = static_cast<T> (zSrcOffset);
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
            size_t ix = it.get_global_id(0);
            size_t iz = it.get_global_id(1);
            auto idst = ::gridToIndex(nGridX, ix, iz);
            T delX = xShiftedSource - ix*hx;
            T delZ = zShiftedSource - iz*hz;
            travelTimesDevice[idst] = sycl::hypot(delX, delZ)*slowness;
        });
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

template<typename T>
void gradient2d(const size_t nGridX, const size_t nGridZ,
                const double xSrcOffset, const double zSrcOffset,
                const double dx, const double dz, 
                const double vel,
                std::vector<T> *gradient)
{
    const T epsilon = std::numeric_limits<T>::epsilon();
    sycl::queue q{sycl::default_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> ();
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    // Set memory on the host
    auto nGrid = nGridX*nGridZ;
    if (gradient->size() != 2*nGrid){gradient->resize(2*nGrid, 0);}
    // Determine ranges
    sycl::range global{nGridX, nGridZ};
    auto nTileX = std::min(workGroupSize, nGridX);
    auto nTileZ = std::min(workGroupSize, nGridZ);
    sycl::range local{nTileX, nTileZ};
    // Simplify some geometric terms
    auto xShiftedSource = static_cast<T> (xSrcOffset);
    auto zShiftedSource = static_cast<T> (zSrcOffset);
    auto hx = static_cast<T> (dx);
    auto hz = static_cast<T> (dz);
    auto slowness = static_cast<T> (1./vel); // Do division here
    // When this goes out of scope the data should be copied back to the host
    {
    // Create the buffers
    sycl::buffer<T> gradientBuffer(*gradient);
    // Tabulate the gradient which involves taking the derivative of a 
    // distance function and scaling by the slowness
    auto eGradient = q.submit([&](sycl::handler &h) 
    {
        sycl::accessor gradientAccessor(gradientBuffer, h,
                                        sycl::write_only, sycl::no_init);
        h.parallel_for(sycl::nd_range{global, local},
                       [=](sycl::nd_item<2> it) 
        {
            size_t ix = it.get_global_id(0);
            size_t iz = it.get_global_id(1);
            auto idst = 2*::gridToIndex(nGridX, ix, iz);
            T delX = ix*hx - xShiftedSource;
            T delZ = iz*hz - zShiftedSource;
            // Let this thing safely go to zero 
            T distance = sycl::fmax(epsilon, sycl::hypot(delX, delZ));
            T invDistanceSlowness = slowness/distance;
            // Analytic expression for gradient
            gradientAccessor[idst]     = delX*invDistanceSlowness;
            gradientAccessor[idst + 1] = delZ*invDistanceSlowness;
        });
    });
    } // Result copied back to host since gradientBuffer goes out of scope
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
    std::vector<T> mTravelTimeGradientField;
    double mVelocity{0};
    bool mHaveTravelTimeField{false};
    bool mHaveTravelTimeGradientField{false};
    bool mHaveSource{false};
    bool mInitialized{false}; 
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

/// Reset the class
template<class T>
void Homogeneous2D<T>::clear() noexcept
{
    pImpl->mGeometry.clear();
    pImpl->mSource.clear();
    pImpl->mTravelTimeField.clear();
    pImpl->mTravelTimeGradientField.clear();
    pImpl->mVelocity = 0;
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveTravelTimeGradientField = false;
    pImpl->mHaveSource = false;
    pImpl->mInitialized = false;
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

/// Initialize the class
template<class T>
void Homogeneous2D<T>::initialize(const Geometry2D &geometry)
{
    clear();
    if (!geometry.haveNumberOfGridPointsInX())
    {   
        throw std::invalid_argument("Grid points in x not set on geometry");
    }   
    if (!geometry.haveNumberOfGridPointsInZ())
    {   
        throw std::invalid_argument("Grid points in z not set on geometry");
    }   
    if (!geometry.haveGridSpacingInX())
    {   
        throw std::invalid_argument("Grid spacing in x not set on geometry");
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
bool Homogeneous2D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Set source
template<class T>
void Homogeneous2D<T>::setSource(const EikonalXX::Source2D &source)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveTravelTimeGradientField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!source.haveLocationInX())
    {   
        throw std::invalid_argument("Source location in x not set");
    }   
    if (!source.haveLocationInZ())
    {   
        throw std::invalid_argument("Source location in z not set");
    }   
    pImpl->mSource = source;
    pImpl->mHaveSource = true;
}

/// Set the source from an (x,z) pair
template<class T>
void Homogeneous2D<T>::setSource(
    const std::pair<double, double> &sourceLocation)
{
    pImpl->mHaveSource = false; 
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveTravelTimeGradientField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
    auto xmin = pImpl->mGeometry.getOriginInX();
    auto zmin = pImpl->mGeometry.getOriginInZ();
    auto xmax = xmin + (nx - 1)*dx;
    auto zmax = zmin + (nz - 1)*dz; 
    if (sourceLocation.first < xmin || sourceLocation.first > xmax)
    {
        throw std::invalid_argument("x source position = "
                                  + std::to_string(sourceLocation.first)
                                  + " must be in range [" + std::to_string(xmin)
                                  + "," + std::to_string(xmax) + "]");
    }
    if (sourceLocation.second < zmin || sourceLocation.second > zmax)
    {
        throw std::invalid_argument("z source position = "
                                  + std::to_string(sourceLocation.second)
                                  + " must be in range [" + std::to_string(zmin)
                                  + "," + std::to_string(zmax) + "]");
    } 
    // Create a source
    EikonalXX::Source2D source;
    source.setGeometry(pImpl->mGeometry);
    source.setLocationInX(sourceLocation.first);
    source.setLocationInZ(sourceLocation.second);
    setSource(source);
}

/// Have source?
template<class T>
bool Homogeneous2D<T>::haveSource() const noexcept
{
    return pImpl->mHaveSource;
}

/// Get source
template<class T>
EikonalXX::Source2D Homogeneous2D<T>::getSource() const
{
    if (!haveSource())
    {
        throw std::runtime_error("Source location not set");
    }
    return pImpl->mSource;
}

/// Set velocity model
template<class T>
void Homogeneous2D<T>::setVelocityModel(const double velocity)
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveTravelTimeGradientField = false;
    if (velocity <= 0)
    {
        throw std::invalid_argument("Velocity = " + std::to_string(velocity)
                                  + " must be positive");
    }
    pImpl->mVelocity = velocity;
}

/// Have velocity?
template<class T>
bool Homogeneous2D<T>::haveVelocityModel() const noexcept
{
    return pImpl->mVelocity > 0;
}

template<class T>
T Homogeneous2D<T>::getSlowness(const int iCellX, const int iCellZ) const 
{
    if (!haveVelocityModel())
    {
         throw std::runtime_error("Velocity model not set");
    }
    auto nCellX = pImpl->mGeometry.getNumberOfCellsInX();
    auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
    if (iCellX < 0 || iCellX >= nCellX)
    {
        throw std::invalid_argument("iCellX = " + std::to_string(iCellX)
                                  + " must be in range [0,"
                                  + std::to_string(nCellX) + "]");
    }
    if (iCellZ < 0 || iCellZ >= nCellZ)
    {
        throw std::invalid_argument("iCellZ = " + std::to_string(iCellZ)
                                  + " must be in range [0,"
                                  + std::to_string(nCellZ) + "]");
    }
    return 1./pImpl->mVelocity;
}

/// Solve the eikonal equation
template<class T>
void Homogeneous2D<T>::solve()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!haveSource()){throw std::runtime_error("Source not yet set");}
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    pImpl->mHaveTravelTimeField = false;
    auto nx = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInX());
    auto nz = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInZ());
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto xSrcOffset = pImpl->mSource.getOffsetInX();
    auto zSrcOffset = pImpl->mSource.getOffsetInZ();
    // Solve it
    ::solve2d(nx, nz,
              xSrcOffset, zSrcOffset,
              dx, dz,
              pImpl->mVelocity,
              &pImpl->mTravelTimeField);
    pImpl->mHaveTravelTimeField = true;
}

template<class T>
void Homogeneous2D<T>::computeTravelTimeGradientField()
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!haveSource()){throw std::runtime_error("Source not yet set");}
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    auto nx = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInX());
    auto nz = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPointsInZ());
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto xSrcOffset = pImpl->mSource.getOffsetInX();
    auto zSrcOffset = pImpl->mSource.getOffsetInZ();
    // Compute gradient
    ::gradient2d(nx, nz,
                 xSrcOffset, zSrcOffset,
                 dx, dz,
                 pImpl->mVelocity,
                 &pImpl->mTravelTimeGradientField);
    pImpl->mHaveTravelTimeGradientField = true; 
}

/// Have travel time field?
template<class T>
bool Homogeneous2D<T>::haveTravelTimeField() const noexcept
{
    return pImpl->mHaveTravelTimeField;
}

template<class T>
std::vector<T> Homogeneous2D<T>::getTravelTimeField() const
{
    if (!haveTravelTimeField())
    {
        throw std::runtime_error("Travel time field not yet computed");
    }
    return pImpl->mTravelTimeField;
}

template<class T>
const T* Homogeneous2D<T>::getTravelTimeFieldPointer() const
{
    if (!haveTravelTimeField())
    {
        throw std::runtime_error("Travel time field not yet computed");
    }
    return pImpl->mTravelTimeField.data();
}

/// Have gradient fields?
template<class T>
bool Homogeneous2D<T>::haveTravelTimeGradientField() const noexcept
{
    return pImpl->mHaveTravelTimeGradientField;
}

template<class T>
std::vector<T> Homogeneous2D<T>::getTravelTimeGradientField() const
{
    if (!haveTravelTimeGradientField())
    {
         throw std::runtime_error("Travel time gradient field not computed");
    }
    return pImpl->mTravelTimeGradientField;
}

template<class T>
const T* Homogeneous2D<T>::getTravelTimeGradientFieldPointer() const
{
    if (!haveTravelTimeGradientField())
    {   
         throw std::runtime_error("Travel time gradient field not computed");
    }   
    return pImpl->mTravelTimeGradientField.data();
}

/// Write the travel time field
template<class T>
void Homogeneous2D<T>::writeVTK(const std::string &fileName,
                                const std::string &title,
                                const bool writeGradient) const
{
    auto tPtr = getTravelTimeFieldPointer(); // Throws
    const T *gradientPtr = nullptr;
    if (writeGradient){gradientPtr = getTravelTimeGradientFieldPointer();}
    IO::VTKRectilinearGrid2D vtkWriter;
    constexpr bool writeBinary{true};
    vtkWriter.open(fileName, pImpl->mGeometry, title, writeBinary); 
    vtkWriter.writeNodalDataset(title, tPtr, EikonalXX::Ordering2D::Natural);
    if (writeGradient)
    {
        vtkWriter.writeNodalVectorDataset(title + "_gradient",
                                          gradientPtr,
                                          EikonalXX::Ordering2D::Natural);
    }
    vtkWriter.close();
}

///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///

template class EikonalXX::Analytic::Homogeneous2D<double>;
template class EikonalXX::Analytic::Homogeneous2D<float>;
