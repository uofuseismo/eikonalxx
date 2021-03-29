#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
#if __has_include(<pstl/execution>)
   #include <pstl/execution>
   #include <pstl/algorithm>
   #define USE_PSTL
#endif
#include "eikonalxx/solver2d.hpp"
#include "private/solver2d.hpp"
#include "private/timer.hpp"
#include "eikonalxx/solverOptions.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/geometry2d.hpp"

using namespace EikonalXX;

template<class T>
class Solver2D<T>::Solver2DImpl
{
public:
    /// Initialize travel times and set travel times around source
    void initializeTravelTimes()
    {
        auto verbosity = mOptions.getVerbosity();
        bool ldebug = (verbosity == Verbosity::DEBUG);
        if (ldebug)
        {
            std::cout << "Initializing travel times near source" << std::endl;
        }
        int nGridX = mGeometry.getNumberOfGridPointsInX();
        int nGridZ = mGeometry.getNumberOfGridPointsInZ();
        const auto sPtr = mVelocityModel.getSlownessPointer();
        auto sourceSlowness = static_cast<T> (sPtr[mSourceCell]);
        auto dx = static_cast<T> (mGeometry.getGridSpacingInX());
        auto dz = static_cast<T> (mGeometry.getGridSpacingInZ());
        auto xShiftedSource = static_cast<T> (mSource.getOffsetInX()); //static_cast<T> (mShiftedSource.first);
        auto zShiftedSource = static_cast<T> (mSource.getOffsetInZ()); //static_cast<T> (mShiftedSource.second);
        std::fill(mTravelTimeField.begin(), mTravelTimeField.end(), mHuge);
        std::vector<std::pair<int, int>> sourceNodes;
        auto sourceCellX = mSource.getCellInX();
        auto sourceCellZ = mSource.getCellInZ();
        sourceNodes.push_back({sourceCellX    , sourceCellZ});
        sourceNodes.push_back({sourceCellX + 1, sourceCellZ});
        sourceNodes.push_back({sourceCellX,     sourceCellZ + 1});
        sourceNodes.push_back({sourceCellX + 1, sourceCellZ + 1});
        for (int i=0; i<static_cast<int> (sourceNodes.size()); ++i)
        {
            auto iSrc = gridToIndex(nGridX,
                                    sourceNodes[i].first,
                                    sourceNodes[i].second);
            mTravelTimeField.at(iSrc)
                = computeAnalyticalTravelTime(sourceNodes[i].first,
                                              sourceNodes[i].second,
                                              dx, dz,
                                              xShiftedSource, zShiftedSource,
                                              sourceSlowness);
            if (ldebug)
            {
                std::cout << "Travel time at node (ix,iz) = ("
                          << sourceNodes[i].first << ","
                          << sourceNodes[i].second << ") is "
                          << mTravelTimeField[iSrc] << " (s)" << std::endl;
            }
        }
        // Initialize update nodes then freeze nodes around source
        mSolverSweep1.initializeUpdateNodes(sourceNodes, verbosity);
        mSolverSweep2.initializeUpdateNodes(sourceNodes, verbosity);
        mSolverSweep3.initializeUpdateNodes(sourceNodes, verbosity);
        mSolverSweep4.initializeUpdateNodes(sourceNodes, verbosity);
    }
    /// The velocity model
    Model2D<T> mVelocityModel;
    /// The model geometry
    Geometry2D mGeometry;
    /// The solver options
    SolverOptions mOptions;
    /// The source
    Source2D mSource;
    /// Solver for each sweep direction
    Solver2DSweep<T, SweepNumber2D::SWEEP1> mSolverSweep1;
    Solver2DSweep<T, SweepNumber2D::SWEEP2> mSolverSweep2; 
    Solver2DSweep<T, SweepNumber2D::SWEEP3> mSolverSweep3;
    Solver2DSweep<T, SweepNumber2D::SWEEP4> mSolverSweep4;
    /// The travel time field.
    std::vector<T> mTravelTimeField;
    /// The source location
    //std::pair<double, double> mSourceLocation{0, 0};
    /// The shifted source location in (x,z).  This is useful to the solver.
    //std::pair<T, T> mShiftedSource{0, 0};
    /// Huge value for travel time initialization
    const T mHuge = HUGE;
    /// The source cell index in x.
    //int mSourceCellX = 0;
    /// The source cell index in z.
    //int mSourceCellZ = 0;
    /// The source cell - this is used to extract the slowness at the source.
    int mSourceCell = 0;
    /// Initialized?
    bool mInitialized = false;
    /// Have velocity model?
    bool mHaveVelocityModel = false;
    /// Have travel time field?
    bool mHaveTravelTimeField = false;
    /// Uniform grid?
    bool mUniformGrid = false;
    /// Have the source?
    bool mHaveSource = false;
};

/// C'tor
template<class T>
Solver2D<T>::Solver2D() : 
    pImpl(std::make_unique<Solver2DImpl> ())
{
}

/// Destructor
template<class T>
Solver2D<T>::~Solver2D() = default;

/// Reset the class
template<class T>
void Solver2D<T>::clear() noexcept
{
    pImpl->mVelocityModel.clear();
    pImpl->mGeometry.clear();
    pImpl->mOptions.clear();
    pImpl->mSource.clear();
    pImpl->mSolverSweep1.clear();
    pImpl->mSolverSweep2.clear();
    pImpl->mSolverSweep3.clear();
    pImpl->mSolverSweep4.clear();
    pImpl->mTravelTimeField.clear();
    //pImpl->mSourceLocation = std::make_pair<double, double> (0, 0);
    //pImpl->mShiftedSource = std::make_pair<T, T> (0, 0);
    //pImpl->mSourceCellX = 0;
    //pImpl->mSourceCellZ = 0;
    pImpl->mSourceCell = 0;
    pImpl->mHaveVelocityModel = false;
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    pImpl->mUniformGrid = false;
    pImpl->mInitialized = false;
}

/// Initializes the solver
template<class T>
void Solver2D<T>::initialize(const SolverOptions &options,
                             const Geometry2D &geometry)
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
    pImpl->mOptions = options;
    // Initialize the solvers
    pImpl->mSolverSweep1.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep2.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep3.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep4.initialize(pImpl->mGeometry, pImpl->mOptions);
    // Allocate space for travel time field
    pImpl->mTravelTimeField.resize(pImpl->mGeometry.getNumberOfGridPoints(), 0);
    pImpl->mInitialized = true;
    // Print some debug info
    if (pImpl->mOptions.getVerbosity() >= Verbosity::INFO)
    {
        std::cout << "Number of grid points in x: "
                  << pImpl->mSolverSweep1.mGridX << std::endl;
        std::cout << "Number of grid points in z: "
                  << pImpl->mSolverSweep1.mGridZ << std::endl;
        if (pImpl->mSolverSweep1.mUniformGrid)
        {
            std::cout << "Uniform grid spacing: " << pImpl->mSolverSweep1.mDx
                      << " (m)" << std::endl;
        }
        else
        {
            std::cout << "Grid spacing in x: " << pImpl->mSolverSweep1.mDx
                      << " (m)" << std::endl;
            std::cout << "Grid spacing in z: " << pImpl->mSolverSweep1.mDz
                      << " (m)" << std::endl;
        }
    }
}

/// Initialized?
template<class T>
bool Solver2D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Gets the geometry
template<class T>
Geometry2D Solver2D<T>::getGeometry() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGeometry;
}

/// Gets the model options
template<class T>
SolverOptions Solver2D<T>::getOptions() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mOptions;
}

template<class T>
void Solver2D<T>::setSource(const Source2D &source)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!source.haveLocationInX())
    {
        throw std::invalid_argument("Source location in x not set");
    }
    if (!source.haveLocationInZ())
    {
        throw std::invalid_argument("Source location in z not set");
    }
    //pImpl->mShiftedSource = std::make_pair<T, T> (source.getOffsetInX(),
    //                                              source.getOffsetInZ());
    //pImpl->mSourceCellX = source.getSourceCellInX();
    //pImpl->mSourceCellZ = source.getSourceCellInZ(); 
    pImpl->mSourceCell = gridToIndex(pImpl->mGeometry.getNumberOfCellsInX(),
                                     source.getCellInX(), source.getCellInZ());
    //                                 pImpl->mSourceCellX, pImpl->mSourceCellZ);
    pImpl->mSource = source;
    pImpl->mHaveSource = true;
    if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
    {
        std::cout << pImpl->mSource << std::endl;
    } 
}

template<class T>
void Solver2D<T>::setSource(const std::pair<double, double> &sourceLocation)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
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
    // Tell user what is about to happen
    if (pImpl->mOptions.getVerbosity() >= Verbosity::INFO)
    {
        std::cout << "Setting source location (x,z)=(" << sourceLocation.first
                  << "," << sourceLocation.second << ")" << std::endl;
    }
    // Get some solver information 
    Source2D source;
    source.setGeometry(pImpl->mGeometry);
    source.setLocationInX(sourceLocation.first);
    source.setLocationInZ(sourceLocation.second);
    setSource(source);
/*
    pImpl->mShiftedSource = std::make_pair<T, T> (sourceLocation.first  - xmin,
                                                  sourceLocation.second - zmin);
    pImpl->mSourceCellX = static_cast<int> (pImpl->mShiftedSource.first/dx);
    pImpl->mSourceCellZ = static_cast<int> (pImpl->mShiftedSource.second/dz);
    pImpl->mSourceCell = gridToIndex(pImpl->mGeometry.getNumberOfCellsInX(),
                                     pImpl->mSourceCellX, pImpl->mSourceCellZ);
    pImpl->mSourceLocation = sourceLocation;
    pImpl->mHaveSource = true;
    // Some fine-grained debug
    if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
    {
        std::cout << "Shifted source location (x,z)=("
                  << pImpl->mShiftedSource.first << "," 
                  << pImpl->mShiftedSource.second << ")" << std::endl;
        std::cout << "Source cell (iCellX,iCellZ)=(" 
                  << pImpl->mSourceCellX << "," << pImpl->mSourceCellZ
                  << ")" << std::endl;
    }
*/
}

/// Get source
template<class T>
Source2D Solver2D<T>::getSource() const
{
    if (!haveSource())
    {
        throw std::runtime_error("Source location not set");
    }
    return pImpl->mSource;
}

/// Have source location?
template<class T>
bool Solver2D<T>::haveSource() const noexcept
{
    return pImpl->mHaveSource;
}

/// Sets the velocity model
template<class T>
void Solver2D<T>::setVelocityModel(const Model2D<T> &velocityModel)
{
    pImpl->mHaveVelocityModel = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    if (!velocityModel.isInitialized())
    {
        throw std::invalid_argument("Velocity model not initialized");
    }
    if (velocityModel.getGeometry() != pImpl->mGeometry)
    {
        throw std::invalid_argument(
            "Velocity model's geometry does not match the solver's geometry");
    }
    if (pImpl->mOptions.getVerbosity() >= Verbosity::INFO)
    {
        std::cout << "Setting velocity model..." << std::endl;
    } 
    // Copy the velocity model
    pImpl->mVelocityModel = velocityModel;
    // Set the velocity models for each sweep
    pImpl->mSolverSweep1.setVelocityModel(velocityModel);
    pImpl->mSolverSweep2.setVelocityModel(velocityModel);
    pImpl->mSolverSweep3.setVelocityModel(velocityModel);
    pImpl->mSolverSweep4.setVelocityModel(velocityModel);
    pImpl->mHaveVelocityModel = true;
}

/// Gets the velocity model
template<class T>
Model2D<T> Solver2D<T>::getVelocityModel() const
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    return pImpl->mVelocityModel;
}

/// Have the velocity model?
template<class T>
bool Solver2D<T>::haveVelocityModel() const noexcept
{
    return pImpl->mHaveVelocityModel;
}

/// Applies the solver
template<class T>
void Solver2D<T>::solve()
{
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (!haveSource())
    {
        throw std::runtime_error("Source location not set");
    }
    // Initialize travel time field based on the source
    Timer timer;
    auto verbosity = pImpl->mOptions.getVerbosity();
    timer.start();
    pImpl->initializeTravelTimes();
    timer.end();
    if (verbosity == Verbosity::DEBUG)
    {
        std::cout << "Travel time field initialization time: "
                  << timer.getDuration() << " (s)" << std::endl;
    }
    // Set source information on solver
    int iSrcX = pImpl->mSource.getCellInX(); //pImpl->mSourceCellX;
    int iSrcZ = pImpl->mSource.getCellInZ(); //pImpl->mSourceCellZ;
    T xSourceOffset = static_cast<T> (pImpl->mSource.getOffsetInX()); //pImpl->mShiftedSource.first;
    T zSourceOffset = static_cast<T> (pImpl->mSource.getOffsetInZ()); //pImpl->mShiftedSource.second;
    auto slownessPtr = pImpl->mVelocityModel.getSlownessPointer();
    T sourceSlowness = slownessPtr[pImpl->mSourceCell];
    pImpl->mSolverSweep1.setSourceInformation(iSrcX, iSrcZ,
                                              xSourceOffset, zSourceOffset,
                                              sourceSlowness);
    pImpl->mSolverSweep2.setSourceInformation(iSrcX, iSrcZ,
                                              xSourceOffset, zSourceOffset,
                                              sourceSlowness);
    pImpl->mSolverSweep3.setSourceInformation(iSrcX, iSrcZ,
                                              xSourceOffset, zSourceOffset,
                                              sourceSlowness);
    pImpl->mSolverSweep4.setSourceInformation(iSrcX, iSrcZ,
                                              xSourceOffset, zSourceOffset,
                                              sourceSlowness);
    // Get number of iterations (refinements) 
    auto nIterations = pImpl->mOptions.getNumberOfSweeps();
    // Perform fast sweeping on CPU
    auto algorithm = pImpl->mOptions.getAlgorithm();
    if (algorithm == SolverAlgorithm::FAST_SWEEPING_METHOD)
    {
        if (verbosity == Verbosity::DEBUG)
        {
            std::cout << "Initializing fast sweeping method..." << std::endl;
        }
        // Initialization sweeps
        timer.start(); 
        auto tTimesPtr = pImpl->mTravelTimeField.data();
        pImpl->mSolverSweep1.initializationFSM(slownessPtr, tTimesPtr);
std::cout << "---------------------------------------------------------------" << std::endl;
        pImpl->mSolverSweep2.initializationFSM(slownessPtr, tTimesPtr);
std::cout << "---------------------------------------------------------------" << std::endl;
        pImpl->mSolverSweep3.initializationFSM(slownessPtr, tTimesPtr);
std::cout << "---------------------------------------------------------------" << std::endl;
        pImpl->mSolverSweep4.initializationFSM(slownessPtr, tTimesPtr);
        if (verbosity == Verbosity::DEBUG)
        {
            std::cout << "Initialization time: "
                      << timer.getDuration() << " (s)" << std::endl; 
        }
        timer.end();
        // Refinement sweeps
        for (int k=0; k<nIterations; ++k)
        {
        }
    }
    else // Perform level set method on device
    {
        if (verbosity == Verbosity::DEBUG)
        {
            std::cout << "Selecting queue..." << std::endl;
        }
        sycl::queue q{sycl::cpu_selector()};
        // Set space
        auto maxNodesInLevel = pImpl->mSolverSweep1.mMaxNodesInLevel;
        auto nWork = padLength(maxNodesInLevel, sizeof(T), 64);
        T *s0 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *s1 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *s3 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *t0 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *t1 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *t2 = sycl::malloc_device<T> (nWork*sizeof(T), q);
        T *t3 = sycl::malloc_device<T> (nWork*sizeof(T), q);

        q.wait();

        for (int k=0; k<nIterations; ++k)
        {
        }
        q.wait();
        // Release memory
        sycl::free(s0, q);
        sycl::free(s1, q);
        sycl::free(s3, q);
        sycl::free(t0, q);
        sycl::free(t1, q);
        sycl::free(t2, q);
        sycl::free(t3, q);
    }
    
/*
    // Create a queue
    // Get information about geometry, sweeps, and source
    int nLevels = pImpl->mLevels; 
    int nSweepDirections = pImpl->mSweepDirections;
    int sphericalRadius = pImpl->mOptions.getSphericalSolverRadius();
    int nSweeps = pImpl->mOptions.getNumberOfSweeps();
    auto nGrid = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPoints());
    int nGridX = pImpl->mGeometry.getNumberOfGridPointsInX();
    int nGridZ = pImpl->mGeometry.getNumberOfGridPointsInZ();
    auto iSourceCellX = pImpl->mSourceCellX;
    auto iSourceCellZ = pImpl->mSourceCellZ;
    auto dx = static_cast<T> (pImpl->mGeometry.getGridSpacingInX());
    auto dz = static_cast<T> (pImpl->mGeometry.getGridSpacingInZ());
    auto huge = pImpl->mHuge;
    auto maxNodesInLevel = pImpl->mMaxNodesInLevel;
    auto xShiftedSource = static_cast<T> (pImpl->mShiftedSource.first);
    auto zShiftedSource = static_cast<T> (pImpl->mShiftedSource.second);
    const auto sPtr = pImpl->mVelocityModel.getSlownessPointer();
    auto sourceSlowness = sPtr[pImpl->mSourceCell];
    const auto levelOffset = pImpl->mLevelOffset.data();
    bool uniformGrid = pImpl->mUniformGrid;
    // Precompute some constants for the solver
    auto dx2 = pImpl->mDx2;
    auto dz2 = pImpl->mDz2;
    auto dx_dz = pImpl->mDxDividedByDz;
    auto dz_dx = pImpl->mDzDividedByDx;
    auto dx2_p_dz2_inv = pImpl->mDx2PlusDz2Inv;
    auto cosTheta = pImpl->mCosTheta;
    auto sinTheta = pImpl->mSinTheta;
    T dtMax = 0;
    // Create unified shared memory to hold slownesses/travel times in sweep
    auto nWork = padLength(maxNodesInLevel, sizeof(T), 64);
    T *s0 = sycl::malloc_device<T> (nWork*sizeof(T), q); 
    T *s1 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    T *s3 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    T *t0 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    T *t1 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    T *t2 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    T *t3 = sycl::malloc_device<T> (nWork*sizeof(T), q);
    int8_t *updateNode = sycl::malloc_device<int8_t> (nWork*sizeof(int8_t), q);
    int converged = 0;
    sycl::buffer convergedBuffer(&converged, sycl::range(1));
    nWork = padLength(nGrid, sizeof(T), 64);
    T *tTimes = sycl::malloc_shared<T> (nWork*sizeof(T), q);
    // Initialize the traveltimes to huge and set the traveltimes around source
    pImpl->initializeTravelTimes(); 
    auto travelTimePtr = pImpl->mTravelTimeField.data();
    // Copy these traveltimes onto shared pointer
    auto eSetTravelTimes = q.memcpy(tTimes, travelTimePtr, nGrid*sizeof(T));
    // Initialization loop
    for (int sweep=0; sweep<nSweepDirections; ++sweep)
    {
        const auto sweepSlownessPtr = pImpl->getSweepSlownessPointer(sweep);
        const auto updateNodePtr = pImpl->getUpdateNodePointer(sweep);
        if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
        {
            timer.start();
            std::cout << "Performing initialization sweep: "
                      << sweep << std::endl; 
        }
        for (int level=0; level<nLevels; ++level)
        {
            auto offset = levelOffset[level];
            int i0, i1;
            getLevelStartStopIndices(nGridX, nGridZ, level, &i0, &i1);
            auto nNodesInLevel = static_cast<size_t> (i1 - i0);
            sycl::range<1> nodesInLevel{nNodesInLevel};
            // Begin migrating slownesses
            auto eCopySlowness0 = q.memcpy(s0, sweepSlownessPtr[level].s0,
                                           nNodesInLevel*sizeof(T));
            auto eCopySlowness1 = q.memcpy(s1, sweepSlownessPtr[level].s1,
                                           nNodesInLevel*sizeof(T));
            auto eCopySlowness3 = q.memcpy(s3, sweepSlownessPtr[level].s3,
                                           nNodesInLevel*sizeof(T));
            auto eCopyUpdateNodes = q.memcpy(updateNode,
                                             updateNodePtr + offset,
                                             nNodesInLevel*sizeof(int8_t));
            // Parallel loop on level
            auto eGetTravelTimes = q.submit([&](sycl::handler &h)
            {
                h.depends_on(eSetTravelTimes);
                h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                {
                    int it0, it1, it2, it3;
                    sweepLevelIndexToTravelTimeIndices(sweep, level, i0 + i,
                                                       nGridX, nGridZ,
                                                       &it0, &it1, &it2, &it3);
                    t0[i] = tTimes[it0];
                    t1[i] = tTimes[it1];
                    t2[i] = tTimes[it2];
                    t3[i] = tTimes[it3];
                });
            });
            // Compute the finite difference updates
            sycl::event eFiniteDifference;
            if (uniformGrid)
            {
                eFiniteDifference = q.submit([&](sycl::handler &h)
                {
                    h.depends_on(eCopySlowness0);
                    h.depends_on(eCopySlowness1);
                    h.depends_on(eCopySlowness3);
                    h.depends_on(eGetTravelTimes);
                    h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                    {
                        int ix, iz;
                        sweepLevelIndexToGrid(sweep, level, i0 + i,
                                              nGridX, nGridZ,
                                              &ix, &iz);
                        T tUpd =  sphericalFiniteDifference(sphericalRadius,
                                                            huge,
                                                            dx,
                                                            sourceSlowness,
                                                            xShiftedSource,
                                                            zShiftedSource,
                                                            sweep,
                                                            ix, iz,
                                                            s0[i], s1[i], s3[i],
                                                            t1[i], t2[i], t3[i]);
                        //T tUpd = cartesianFiniteDifference(huge,
                        //                                   dx,
                        //                                   s0[i], s1[i], s3[i],
                        //                                   t1[i], t2[i], t3[i]);
                        t0[i] = sycl::fmin(t0[i], tUpd);
                    });
                });
            }
            else
            {
                eFiniteDifference = q.submit([&](sycl::handler &h)
                {
                    h.depends_on(eCopySlowness0);
                    h.depends_on(eCopySlowness1);
                    h.depends_on(eCopySlowness3);
                    h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                    {
                        int ix, iz;
                        sweepLevelIndexToGrid(sweep, level, i0 + i,
                                              nGridX, nGridZ,
                                              &ix, &iz);
                        T tUpd = cartesianFiniteDifference(huge,
                                                           dx, dz,
                                                           dx_dz, dz_dx,
                                                           cosTheta, sinTheta,
                                                           s0[i], s1[i], s3[i],
                                                           t1[i], t2[i], t3[i]);
                        t0[i] = sycl::fmin(t0[i], tUpd);
                    });
                });
            }
            // Copy the updated travel times back to the host
            q.submit([&](sycl::handler &h)
            {
                h.depends_on(eFiniteDifference);
                h.depends_on(eCopyUpdateNodes);
                h.parallel_for(nodesInLevel,
                               [=](sycl::id<1> i)
                {
                    auto j = sweepLevelIndexToIndex(sweep, level, i0 + i,
                                                    nGridX, nGridZ);
                    if (updateNode[i] == UPDATE_NODE){tTimes[j] = t0[i];}
                });
            });
            // Block until I'm done updating travel times in this level
            q.wait();
        } // Loop on levels
    }
    q.wait_and_throw(); 
    // Perform Gauss-Seidel iterations 
    for (int k=0; k<nSweeps; ++k)
    {
        if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
        {
            timer.start();
            std::cout << "Beginning iteration: " << k + 1 << std::endl;
        }
        for (int sweep=0; sweep<nSweepDirections; ++sweep)
        {
            const auto sweepSlownessPtr = pImpl->getSweepSlownessPointer(sweep);
            const auto updateNodePtr = pImpl->getUpdateNodePointer(sweep);
            // Loop on levels in sweep
            for (int level=0; level<nLevels; ++level)
            {
                auto offset = levelOffset[level];
                int i0, i1; 
                getLevelStartStopIndices(nGridX, nGridZ, level, &i0, &i1);
                auto nNodesInLevel = static_cast<size_t> (i1 - i0);
                sycl::range<1> nodesInLevel{nNodesInLevel};
                // Begin migrating slownesses
                auto eCopySlowness0 = q.memcpy(s0, sweepSlownessPtr[level].s0,
                                               nNodesInLevel*sizeof(T));
                auto eCopySlowness1 = q.memcpy(s1, sweepSlownessPtr[level].s1,
                                               nNodesInLevel*sizeof(T));
                auto eCopySlowness3 = q.memcpy(s3, sweepSlownessPtr[level].s3,
                                               nNodesInLevel*sizeof(T));
                auto eCopyUpdateNodes = q.memcpy(updateNode,
                                                 updateNodePtr + offset,
                                                 nNodesInLevel*sizeof(int8_t));
                // Parallel loop on level
                auto eGetTravelTimes = q.submit([&](sycl::handler &h)
                {
                    h.depends_on(eSetTravelTimes);
                    h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                    {
                        int it0, it1, it2, it3;
                        sweepLevelIndexToTravelTimeIndices(sweep, level, i0 + i,
                                                        nGridX, nGridZ,
                                                        &it0, &it1, &it2, &it3);
                        t0[i] = tTimes[it0];
                        t1[i] = tTimes[it1];
                        t2[i] = tTimes[it2];
                        t3[i] = tTimes[it3]; 
                    }); 
                }); 
                // Compute the finite difference updates
                sycl::event eFiniteDifference;
                if (uniformGrid)
                {
                    eFiniteDifference = q.submit([&](sycl::handler &h)
                    {
                        h.depends_on(eCopySlowness0);
                        h.depends_on(eCopySlowness1);
                        h.depends_on(eCopySlowness3);
                        h.depends_on(eGetTravelTimes);
                        h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                        {
                            T tUpd = cartesianFiniteDifference(huge,
                                                           dx, 
                                                           s0[i], s1[i], s3[i],
                                                           t1[i], t2[i], t3[i]);
                            t0[i] = sycl::fmin(t0[i], tUpd);
                        });
                    });
                }
                else
                {
                    eFiniteDifference = q.submit([&](sycl::handler &h)
                    {
                        h.depends_on(eCopySlowness0);
                        h.depends_on(eCopySlowness1);
                        h.depends_on(eCopySlowness3);
                        h.parallel_for(nodesInLevel, [=](sycl::id<1> i)
                        {
                            T tUpd = cartesianFiniteDifference(huge,
                                                           dx, dz,
                                                           dx_dz, dz_dx,
                                                           cosTheta, sinTheta,
                                                           s0[i], s1[i], s3[i],
                                                           t1[i], t2[i], t3[i]);
                            t0[i] = sycl::fmin(t0[i], tUpd);
                        });
                    }); 
                }
                // Copy the updated travel times back to the host
                q.submit([&](sycl::handler &h)
                {
                    h.depends_on(eFiniteDifference);
                    h.depends_on(eCopyUpdateNodes);
                    h.parallel_for(nodesInLevel,
                                   [=](sycl::id<1> i)
                    {
                        auto j = sweepLevelIndexToIndex(sweep, level, i0 + i,
                                                        nGridX, nGridZ);
                        if (updateNode[i] == UPDATE_NODE){tTimes[j] = t0[i];}
                    });
                });
                // Block until I'm done updating travel times in this level
                q.wait();
            } // Loop on levels
        } // Loop on sweep directions
        if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
        {
            timer.end();
            std::cout << "Iteration took: " << timer.getDuration()
                      << " (s)" << std::endl;
        }
    } // Loop on iterations
    // Copy result back
    q.memcpy(travelTimePtr, tTimes, nGrid*sizeof(T)); 
    q.wait();
std::cout << travelTimePtr[0] << std::endl;
std::cout<< travelTimePtr[50] << std::endl;
std::cout << travelTimePtr[(nGridX-1)*(nGridZ-1)] << std::endl;
    // Release memory
    sycl::free(s0, q);
    sycl::free(s1, q);
    sycl::free(s3, q);
    sycl::free(t0, q);
    sycl::free(t1, q);
    sycl::free(t2, q);
    sycl::free(t3, q);
    sycl::free(updateNode, q);
    sycl::free(tTimes, q);
    pImpl->mHaveTravelTimeField = true;
*/
}

/// Get travel time field
template<class T>
std::vector<T> Solver2D<T>::getTravelTimeField() const
{
    return pImpl->mTravelTimeField;;
}

/// Get travel time field pointer
template<class T>
const T* Solver2D<T>::getTravelTimeFieldPointer() const
{
    if (!haveTravelTimeField())
    {
        throw std::runtime_error("Travel time field not yet computed");
    }
    return pImpl->mTravelTimeField.data();
}

/// Have travel time field?
template<class T>
bool Solver2D<T>::haveTravelTimeField() const noexcept
{
    return pImpl->mHaveTravelTimeField;
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template class EikonalXX::Solver2D<double>;
template class EikonalXX::Solver2D<float>;
