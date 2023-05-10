#include <vector>
#include <limits>
#include <string>
#include <algorithm>
//#include <spdlog/spdlog.h>
//#include <spdlog/sinks/stdout_sinks.h>
//#include <spdlog/sinks/stdout_color_sinks.h>
#ifndef NDEBUG
  #include <cassert>
#endif
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
    Solver2DImpl() //:
        //mStdOutSink(std::make_shared<spdlog::sinks::stdout_color_sink_mt> ())//,
        //mLogger(std::make_shared<spdlog::logger> ("stdout", mStdOutSink))
    {
        //mLogger->set_level(spdlog::level::err);
    } 
    /// Sets the log level
    void setLogLevel(const Verbosity verbosity)
    {
        if (verbosity == Verbosity::Debug)
        {
            //mLogger->set_level(spdlog::level::debug);
        }
        else if (verbosity == Verbosity::Info)
        {
            //mLogger->set_level(spdlog::level::info);
        } 
        else if (verbosity == Verbosity::Warning)
        {
            //mLogger->set_level(spdlog::level::warn);
        }
        else if (verbosity == Verbosity::Error)
        {
            //mLogger->set_level(spdlog::level::err);
        }
    }
    /// Initialize travel times and set travel times around source
    void initializeTravelTimes()
    {
        auto verbosity = mOptions.getVerbosity();
        bool ldebug = (verbosity == Verbosity::Debug);
        //mLogger->debug("Initializing travel times near source.");
        int nGridX = mGeometry.getNumberOfGridPointsInX();
        const auto sPtr = mVelocityModel.getSlownessPointer();
        auto sourceSlowness = static_cast<T> (sPtr[mSourceCell]);
        auto dx = static_cast<T> (mGeometry.getGridSpacingInX());
        auto dz = static_cast<T> (mGeometry.getGridSpacingInZ());
        auto xShiftedSource = static_cast<T> (mSource.getOffsetInX());
        auto zShiftedSource = static_cast<T> (mSource.getOffsetInZ());
        std::fill(mTravelTimeField.begin(), mTravelTimeField.end(), HUGE);
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
    /// The logger
    //std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> mStdOutSink{nullptr};
    //std::shared_ptr<spdlog::logger> mLogger{nullptr};
    /// The velocity model
    Model2D<T> mVelocityModel;
    /// The model geometry
    Geometry2D mGeometry;
    /// The solver options
    SolverOptions mOptions;
    /// The source
    Source2D mSource;
    /// Solver for each sweep direction
    Solver2DSweep<T, SweepNumber2D::Sweep1> mSolverSweep1;
    Solver2DSweep<T, SweepNumber2D::Sweep2> mSolverSweep2; 
    Solver2DSweep<T, SweepNumber2D::Sweep3> mSolverSweep3;
    Solver2DSweep<T, SweepNumber2D::Sweep4> mSolverSweep4;
    /// The travel time field.
    std::vector<T> mTravelTimeField;
    /// The source location
    //std::pair<double, double> mSourceLocation{0, 0};
    /// The shifted source location in (x,z).  This is useful to the solver.
    //std::pair<T, T> mShiftedSource{0, 0};
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
    pImpl->setLogLevel(options.getVerbosity());
    // Initialize the solvers
    pImpl->mSolverSweep1.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep2.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep3.initialize(pImpl->mGeometry, pImpl->mOptions);
    pImpl->mSolverSweep4.initialize(pImpl->mGeometry, pImpl->mOptions);
    // Allocate space for travel time field
    pImpl->mTravelTimeField.resize(pImpl->mGeometry.getNumberOfGridPoints(), 0);
    pImpl->mInitialized = true;
    // Print some debug info
    //pImpl->mLogger->debug("Number of grid points in x: "
    //                    + std::to_string(pImpl->mSolverSweep1.mGridX)); 
    //pImpl->mLogger->debug("Number of grid points in z: "
    //                    + std::to_string(pImpl->mSolverSweep1.mGridZ));
    if (pImpl->mSolverSweep1.mUniformGrid)
    {
        //pImpl->mLogger->debug("Uniform grid spacing of: "
        //                     + std::to_string(pImpl->mSolverSweep1.mDx)
        //                    + " (m)");
    }
    else
    {
        //pImpl->mLogger->debug("Grid spacing in x: "
        //                    + std::to_string(pImpl->mSolverSweep1.mDx)
        //                    + " (m)"); 
        //pImpl->mLogger->debug("Grid spacing in z: "
        //                    + std::to_string(pImpl->mSolverSweep1.mDx)
        //                    + " (m)"); 
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
    pImpl->mSourceCell = source.getCell();
    pImpl->mSource = source;
    pImpl->mHaveSource = true;
    if (pImpl->mOptions.getVerbosity() == Verbosity::Debug)
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
    if (pImpl->mOptions.getVerbosity() >= Verbosity::Info)
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
    //pImpl->mLogger->debug("Setting velocity model...");
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
    timer.start();
    pImpl->initializeTravelTimes();
    timer.end();
    //pImpl->mLogger->debug("Travel time field initialization time: "
    //                    + std::to_string(timer.getDuration()) + " (s)");
    // Set source information on solver
    int iSrcX = pImpl->mSource.getCellInX();
    int iSrcZ = pImpl->mSource.getCellInZ();
    T xSourceOffset = static_cast<T> (pImpl->mSource.getOffsetInX());
    T zSourceOffset = static_cast<T> (pImpl->mSource.getOffsetInZ());
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
    constexpr bool initialize = true;
    constexpr bool noInitialize = false;
    double initializationTime = 0;
    // Perform initialization sweeps
    T *__restrict__ tTimesPtr = pImpl->mTravelTimeField.data();
    //pImpl->mLogger->debug("Initializing fast sweeping method...");
    timer.start(); 
    pImpl->mSolverSweep1.updateFSM(slownessPtr, tTimesPtr, initialize);
//std::cout << "---------------------------------------------------------------" << std::endl;
//getchar();
    pImpl->mSolverSweep2.updateFSM(slownessPtr, tTimesPtr, initialize);
//std::cout << "---------------------------------------------------------------" << std::endl;
    pImpl->mSolverSweep3.updateFSM(slownessPtr, tTimesPtr, initialize);
//std::cout << "---------------------------------------------------------------" << std::endl;
    pImpl->mSolverSweep4.updateFSM(slownessPtr, tTimesPtr, initialize);
    //pImpl->mLogger->debug("Initialization time: "
    //                    + std::to_string(timer.getDuration()) + " (s)");
    timer.end();
    initializationTime = timer.getDuration();

    // Continue with fast sweeping method
    if (algorithm == SolverAlgorithm::FastSweepingMethod)
    {
        // Refinement sweeps
        timer.start();
        for (int k = 0; k < nIterations; ++k)
        {
            pImpl->mSolverSweep1.updateFSM(slownessPtr, tTimesPtr, noInitialize);
            pImpl->mSolverSweep2.updateFSM(slownessPtr, tTimesPtr, noInitialize);
            pImpl->mSolverSweep3.updateFSM(slownessPtr, tTimesPtr, noInitialize);
            pImpl->mSolverSweep4.updateFSM(slownessPtr, tTimesPtr, noInitialize);
        }
    }
    else // Perform level set method on device
    {
        //pImpl->mLogger->debug("Selecting queue...");
        sycl::queue q{sycl::cpu_selector_v};
        q.wait();
        // Refinment sweeps
        timer.start(); 
        for (int k = 0; k < nIterations; ++k)
        {
            pImpl->mSolverSweep1.updateLSM(tTimesPtr, q);
            pImpl->mSolverSweep2.updateLSM(tTimesPtr, q);
            pImpl->mSolverSweep3.updateLSM(tTimesPtr, q);
            pImpl->mSolverSweep4.updateLSM(tTimesPtr, q);
        }
        q.wait();
    }
    timer.end();
    auto updateDuration = timer.getDuration();
    //pImpl->mLogger->debug("Update time: "
    //                    + std::to_string(updateDuration) + " (s)");
    //pImpl->mLogger->debug("Total time: "
    //                    + std::to_string(updateDuration
    //                                   + initializationTime) + " (s)");
    pImpl->mHaveTravelTimeField = true;

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
        if (pImpl->mOptions.getVerbosity() == Verbosity::Debug)
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
        if (pImpl->mOptions.getVerbosity() == Verbosity::Debug)
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
        if (pImpl->mOptions.getVerbosity() == Verbosity::Debug)
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
    return pImpl->mTravelTimeField;
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
