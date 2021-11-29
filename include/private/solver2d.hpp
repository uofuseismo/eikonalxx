#ifndef EIKONALXX_SOLVER2D
#define EIKONALXX_SOLVER2D
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <CL/sycl.hpp>
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include "eikonalxx/enums.hpp"
#include "solverUtilities2d.hpp"

#include <iomanip>

#define HUGE 1.e10
namespace
{

/// Solver class
template<class T, EikonalXX::SweepNumber2D E>
class Solver2DSweep
{
public:
    /// Initialize the variables in this class 
    void initialize(const EikonalXX::Geometry2D &geometry,
                    const EikonalXX::SolverOptions &options)
    {
        // Isotropic grid?
        auto dx = geometry.getGridSpacingInX();
        auto dz = geometry.getGridSpacingInZ();
        if (std::abs(dx - dz) < std::numeric_limits<T>::epsilon()*100)
        {
            mUniformGrid = true;
        }
        // Solver parameters
        mAlgorithm = options.getAlgorithm();
        mFactoredEikonalSolverRadius
            = options.getFactoredEikonalEquationSolverRadius();
        // Geometric information
        mDx = static_cast<T> (dx);
        mDz = static_cast<T> (dz);
        mDxDividedByDz = static_cast<T> (dx/dz);
        mDzDividedByDx = static_cast<T> (dz/dx);
        mCosTheta = static_cast<T> (dx/std::hypot(dx, dz));
        mSinTheta = static_cast<T> (dz/std::hypot(dx, dz));
        // Model dimensions
        mGridX = geometry.getNumberOfGridPointsInX();
        mGridZ = geometry.getNumberOfGridPointsInZ();
        mGrid  = geometry.getNumberOfGridPoints();
        mCellX = geometry.getNumberOfCellsInX();
        mCellZ = geometry.getNumberOfCellsInZ();
        mCell  = geometry.getNumberOfCells();
        // Pointers for LSM solver
        mUpdateNode.resize(mGrid, UPDATE_NODE);
        mLevels = computeNumberOfLevels(mGridX, mGridZ);
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD)
        {
            mLevelOffset = makeLevelOffset(mGridX, mGridZ);
            // Set space for sweep's slowness
            mSweepSlowness.resize(mLevels);
            for (int level=0; level<mLevels; ++level)
            {
                auto nNodes = mLevelOffset[level+1] - mLevelOffset[level];
                mMaxNodesInLevel = std::max(mMaxNodesInLevel, nNodes);
                mSweepSlowness[level].allocate(nNodes);
            }
       }
    }
    /// Initialize the update nodes
    void initializeUpdateNodes(
        const std::vector<std::pair<int, int>> &sourceNodes,
        const EikonalXX::Verbosity verbosity)
    {
        auto ldebug = (verbosity == EikonalXX::Verbosity::DEBUG);
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD)
        {
            auto updateNodePtr = mUpdateNode.data();
            setPreliminaryUpdateNodes<E>(mLevels,
                                         mGridX, mGridZ,
                                         mLevelOffset.data(),
                                         updateNodePtr);
            //auto updatePointer = getUpdateNodePointer(sweep);
            int i0, i1, indx, level, offset;
            for (const auto &sourceNode : sourceNodes)
            {
                gridSweepToLevelIndex<E>(sourceNode.first,
                                         sourceNode.second,
                                         mGridX, mGridZ,
                                         &level, &indx);
                getLevelStartStopIndices(mGridX, mGridZ, level, &i0, &i1);
                offset = mLevelOffset[level];
                if (ldebug)
                {
                    std::cout << "Freezing (sweep,level,index) = ("
                              << static_cast<int> (E) + 1
                              << "," << level << "," << indx
                              << ")" << std::endl;
                }
                auto iDst = offset + (indx - i0);
#ifndef NDEBUG
                assert(iDst >= 0 && iDst < mGridX*mGridZ);
#endif
                updateNodePtr[iDst] = SOURCE_NODE;
            }
        }
        else
        {
            setPreliminaryUpdateNodes<E>(mGridX, mGridZ,
                                         mUpdateNode.data());
            for (const auto &sourceNode : sourceNodes)
            {
                auto indx = gridToIndex(mGridX,
                                        sourceNode.first, sourceNode.second);
                if (ldebug)
                {
                    std::cout << "Freezing (sweep,ix,iz,index) = ("
                              << static_cast<int> (E) + 1 << ","
                              << sourceNode.first << ","
                              << sourceNode.second << ","
                              <<  indx << ")" << std::endl;
                }
                mUpdateNode.at(indx) = SOURCE_NODE;
            }
        }
    }
    /// Clear the class
    void clear() noexcept
    {
        mSweepSlowness.clear();
        mLevelOffset.clear();
        mUpdateNode.clear();
        mDx = 0;
        mDz = 0;
        mDxDividedByDz = 0;
        mDzDividedByDx = 0;
        mCosTheta = 0;
        mSinTheta = 0;
        mSourceOffsetX = 0;
        mSourceOffsetZ = 0;
        mSourceSlowness = 1;
        mLevels = 0;
        mMaxNodesInLevel = 0;
        mGridX = 0;
        mGridZ = 0;
        mGrid = 0;
        mCellX = 0;
        mCellZ = 0;
        mCell = 0;
        mSourceIndexX = 0;
        mSourceIndexZ = 0;
        mFactoredEikonalSolverRadius = 0;
        mAlgorithm = EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD;
        mUniformGrid = false;
    }
    /// Set the velocity model
    void setVelocityModel(const EikonalXX::Model2D<T> &velocityModel)
    {
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD)
        {
#ifndef NDEBUG
            for (int level=0; level<mLevels; ++level)
            {
                mSweepSlowness[level].zero();
            }
#endif
            // Fill slownesses in each sweep
            slownessToSweepSlowness<T, E>(static_cast<size_t> (mLevels),
                                          mGridX, mGridZ, mCellX, mCellZ,
                                          velocityModel.getSlownessPointer(), 
                                          mSweepSlowness.data());
#ifndef NDEBUG
            for (int level=0; level<mLevels; ++level)
            {
                assert(mSweepSlowness[level].getMinimumValue() > 0);
            }
#endif
        }
    }
    /// Sets the source information
    /// @param[in] iSrcX  The source cell in x.
    /// @param[in] iSrcZ  The source cell in z.
    /// @param[in] xSourceOffset   The x source position in meters from the
    ///                            x origin.
    /// @param[in] zSourceOffset   The z source position in meters from the
    ///                            z origin.
    /// @param[in] sourceSlowness  The slowness in s/m in the source cell. 
    void setSourceInformation(const int iSrcX, const int iSrcZ,
                              const T xSourceOffset, const T zSourceOffset,
                              const T sourceSlowness)
    {
        mSourceIndexX = iSrcX;
        mSourceIndexZ = iSrcZ;
        mSourceOffsetX = xSourceOffset;
        mSourceOffsetZ = zSourceOffset;
        mSourceSlowness = sourceSlowness;
    }
    /// @brief Performs the initialization sweep on a grid
    /// @param[in] slowness         The cell-based slowness field in s/m.
    /// @param[in,out] travelTimes  The travel times in at each node in seconds.
    /// @param[in] initialize       True indicates that this is an
    ///                             initialization sweep.
    void updateFSM(const T *__restrict__ slowness,
                   T *__restrict__ travelTimes,
                   const bool initialize)
    {
std::cout.precision(10);
        int ix0, ix1, ixDir, iz0, iz1, izDir;
        if (initialize)
        {
            getLoopLimits<E>(mSourceIndexX, mSourceIndexZ,
                             mGridX, mGridZ,
                             &ix0, &iz0,
                             &ix1, &iz1, 
                             &ixDir, &izDir);
        }
        else
        {
            getLoopLimits<E>(mGridX, mGridZ,
                             &ix0, &iz0,
                             &ix1, &iz1,
                             &ixDir, &izDir);
        }
        int ixShift, izShift, signX, signZ;
        getSweepFiniteDifferenceSigns<E>(&ixShift, &izShift,
                                         &signX, &signZ);
        T t0, t1, t2, t3, tUpd;
        T s0, s1, s3; 
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        T huge = HUGE;
        T diffMax = HUGE;
        if (mUniformGrid)
        {
            auto h = mDx; // dx = dz
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness<E>(ix, iz,
                                                 mCellX, mCellZ,
                                                 &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes<E>(ix, iz, mGridX,
                                                    &it0, &it1, &it2, &it3);
                    t0 = travelTimes[it0];
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    tUpd = finiteDifference(mFactoredEikonalSolverRadius,
                                            huge,
                                            h, 
                                            mSourceSlowness,
                                            ix, iz,
                                            signX, signZ,
                                            ixShift, izShift, 
                                            mSourceIndexX, mSourceIndexZ,
                                            mSourceOffsetX, mSourceOffsetZ,
                                            s0, s1, s3,
                                            t1, t2, t3);
std::cout << ix<< " " << iz << " " << travelTimes[it0] << " " << std::min(travelTimes[it0], tUpd) << std::endl;
                    // Update?
                    //if (mUpdateNode[it0] == UPDATE_NODE)
                    if (tUpd < t0)
                    {
                        travelTimes[it0] = tUpd;
                        diffMax = sycl::fmax(t0 - tUpd, diffMax);
                    }
                }
            }
        }
        else
        {
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                {
                    gridToSurroundingSlowness<E>(ix, iz, 
                                                 mCellX, mCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes<E>(ix, iz, mGridX,
                                                    &it0, &it1, &it2, &it3);
                    t0 = travelTimes[it0];
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    tUpd = finiteDifference(mFactoredEikonalSolverRadius,
                                            huge,
                                            mDx, mDz,
                                            mDxDividedByDz, mDzDividedByDx,
                                            mCosTheta, mSinTheta,
                                            mSourceSlowness,
                                            ix, iz,
                                            signX, signZ,
                                            ixShift, izShift, 
                                            mSourceIndexX, mSourceIndexZ,
                                            mSourceOffsetX, mSourceOffsetZ,
                                            s0, s1, s3, 
                                            t1, t2, t3);
std::cout << "nonuniform: " << ix<< " " << iz << " " << travelTimes[it0] << " " << std::min(travelTimes[it0], tUpd) << std::endl;
                    // Update?
                    //if (mUpdateNode[it0] == UPDATE_NODE)
                    if (tUpd < t0)
                    {
                        travelTimes[it0] = tUpd;
                        diffMax = sycl::fmax(t0 - tUpd, diffMax);
                    }
                }
            }
        } // End check on uniform grid
    }

/*
    /// Performs the fast sweeping method on a grid
    void updateFSM(const T *slowness, T *travelTimes)
    {
        int ix0, ix1, ixDir, iz0, iz1, izDir;
        getLoopLimits(E,
                      mGridX, mGridZ,
                      &ix0, &iz0,
                      &ix1, &iz1, 
                      &ixDir, &izDir);
        int ixShift, izShift, signX, signZ;
        getSweepFiniteDifferenceSigns(E,
                                      &ixShift, &izShift,
                                      &signX, &signZ);
        T t1, t2, t3, tUpd;
        T s0, s1, s3;
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        T huge = HUGE;
        if (mUniformGrid)
        {
            auto h = mDx; // dx = dz
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness(E,
                                              ix, iz,
                                              mCellX, mCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes(E, ix, iz, mGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    ///tUpd = cartesianFiniteDifference(mHuge, h,
                    ///                                 s0, s1, s3,
                    ///                                 t1, t2, t3);
                    tUpd = finiteDifference(mFactoredEikonalSolverRadius,
                                            huge,
                                            h,
                                            mSourceSlowness,
                                            ix, iz,
                                            signX, signZ,
                                            ixShift, izShift,
                                            mSourceIndexX, mSourceIndexZ,
                                            mSourceOffsetX, mSourceOffsetZ,
                                            s0, s1, s3,
                                            t1, t2, t3);
                    // Update?
                    if (mUpdateNode[it0] == UPDATE_NODE)
                    {
                        travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                    }
                }
            }
        }
        else
        {
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                {
                    gridToSurroundingSlowness(E,
                                              ix, iz,
                                              mCellX, mCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes(E, ix, iz, mGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    ///tUpd = cartesianFiniteDifference(huge,
                    ///                                 mDx, mDz,
                    ///                                 mDxDividedByDz,
                    ///                                 mDzDividedByDx,
                    ///                                 mCosTheta, mSinTheta,
                    ///                                 s0, s1, s3,
                    ///                                 t1, t2, t3);

                    tUpd = finiteDifference(mFactoredEikonalSolverRadius,
                                            huge,
                                            mDx, mDz,
                                            mDxDividedByDz, mDzDividedByDx,
                                            mCosTheta, mSinTheta,
                                            mSourceSlowness,
                                            ix, iz, 
                                            signX, signZ,
                                            ixShift, izShift, 
                                            mSourceIndexX, mSourceIndexZ,
                                            mSourceOffsetX, mSourceOffsetZ,
                                            s0, s1, s3, 
                                            t1, t2, t3);
                    // Update?
                    if (mUpdateNode[it0] == UPDATE_NODE)
                    {
                        travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                    }
                }
            }
        } // End check on uniform grid
    } 
*/
//private:
    /// Holds the slowness in the sweep
    std::vector<SweepSlowness2D<T>> mSweepSlowness;
    /// Maps from level to offset in update node
    std::vector<int> mLevelOffset;
    /// Update node for this sweep?
    std::vector<int8_t> mUpdateNode;
    /// Grid spacing in x and z in meters
    T mDx = 0;
    T mDz = 0;
    /// dx/dz
    T mDxDividedByDz = 0;
    /// dz/dx
    T mDzDividedByDx = 0;
    /// dx/sqrt(dx^2 + dz^2)
    T mCosTheta = 0;
    /// dz/sqrt(dx^2 + dz^2)
    T mSinTheta = 0;
    /// x source offset from origin in meters
    T mSourceOffsetX = 0;
    /// z source offset from origin in meters
    T mSourceOffsetZ = 0;
    /// Source slowness (s/m)
    T mSourceSlowness = 1;
    /// Number of levels in level set method
    int mLevels = 0;
    /// Max number of nodes in a level
    int mMaxNodesInLevel = 0;
    /// Number of grid points in x, z, and total
    int mGridX = 0;
    int mGridZ = 0;
    int mGrid = 0;
    /// Number of cells in x, z, and total
    int mCellX = 0;
    int mCellZ = 0;
    int mCell = 0;
    /// Source indices in x and z
    int mSourceIndexX = 0;
    int mSourceIndexZ = 0;
    /// Factored eikonal to cartesian transition
    int mFactoredEikonalSolverRadius = 0;
    /// Algorithm type
    EikonalXX::SolverAlgorithm mAlgorithm
        = EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD;
    /// Is this a uniform grid?
    bool mUniformGrid = false;
};
/*
//----------------------------------------------------------------------------//
//                              Solves +x and +z                              //
//----------------------------------------------------------------------------//
template<class T>
class Solver2DSweep0 //: public Solver2DBaseClass<T>
{
public:
    /// Container with variables for use by solver
    Solver2DVariables<T> mVariables;
    /// Update travel time field with fast sweeping method
    void updateFSM(const T *slowness, T *travelTimes)
    {
        constexpr int sweep = 0;
        //T s0, s1, s3, t1, t2, t3, tUpd;
        auto huge = mVariables.mHuge;
        auto nGridX = mVariables.mGridX;
        auto nGridZ = mVariables.mGridZ;
        auto nCellX = mVariables.mCellX;
        auto nCellZ = mVariables.mCellZ;
        T t1, t2, t3, tUpd;
        T s0, s1, s3;
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        if (mVariables.mUniformGrid)
        {
            auto h = mVariables.mDx;
            for (int iz = 1; iz < nGridZ; ++iz)
            {
                for (int ix = 1; ix < nGridX; ++ix)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    tUpd = cartesianFiniteDifference(huge, h,
                                                     s0, s1, s3,
                                                     t1, t2, t3);
                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        }
        else
        {
            auto dx = mVariables.mDx;
            auto dz = mVariables.mDz;
            auto dx_dz = mVariables.mDxDividedByDz;
            auto dz_dx = mVariables.mDzDividedByDx;
            auto cosTheta = mVariables.mCosTheta;
            auto sinTheta = mVariables.mSinTheta;
            for (int iz = 1; iz < nGridZ; ++iz)
            {
                for (int ix = 1; ix < nGridX; ++ix)
                {
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    tUpd = cartesianFiniteDifference(huge,
                                                     dx, dz,
                                                     dx_dz, dz_dx,
                                                     cosTheta, sinTheta,
                                                     s0, s1, s3,
                                                     t1, t2, t3);

                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        } // End check on uniform grid
    }
};

//----------------------------------------------------------------------------//
//                              Solves -x and +z                              //
//----------------------------------------------------------------------------//
template<class T>
class Solver2DSweep1 //: public Solver2DBaseClass<T>
{
public:
    Solver2DVariables<T> mVariables;
    /// Update travel time field with fast sweeping method
    void updateFSM(const T *slowness, T *travelTimes)
    {
        constexpr int sweep = 1;
        //T s0, s1, s3, t1, t2, t3, tUpd;
        auto huge = mVariables.mHuge;
        auto nGridX = mVariables.mGridX;
        auto nGridZ = mVariables.mGridZ;
        auto nCellX = mVariables.mCellX;
        auto nCellZ = mVariables.mCellZ;
        T t1, t2, t3, tUpd;
        T s0, s1, s3;
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        if (mVariables.mUniformGrid)
        {
            auto h = mVariables.mDx;
            for (int iz = 1; iz < nGridZ; ++iz)
            {
                for (int ix = nGridX - 2; ix >= 0; --ix)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    tUpd = cartesianFiniteDifference(huge, h,
                                                     s0, s1, s3,
                                                     t1, t2, t3);
                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        }
        else
        {
            auto dx = mVariables.mDx;
            auto dz = mVariables.mDz;
            auto dx_dz = mVariables.mDxDividedByDz;
            auto dz_dx = mVariables.mDzDividedByDx;
            auto cosTheta = mVariables.mCosTheta;
            auto sinTheta = mVariables.mSinTheta;
            for (int iz = 1; iz < nGridZ; ++iz)
            {
                for (int ix = nGridX - 2; ix >= 0; --ix)
                {
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    tUpd = cartesianFiniteDifference(huge,
                                                     dx, dz,
                                                     dx_dz, dz_dx,
                                                     cosTheta, sinTheta,
                                                     s0, s1, s3,
                                                     t1, t2, t3);

                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        } // End check on uniform grid
    }
};

//----------------------------------------------------------------------------//
//                              Solves +x and -z                              //
//----------------------------------------------------------------------------//
template<class T>
class Solver2DSweep2 //: public Solver2DBaseClass<T>
{
public:
    Solver2DVariables<T> mVariables;
    /// Update travel time field with fast sweeping method
    void updateFSM(const T *slowness, T *travelTimes)
    {
        constexpr int sweep = 2;
        //T s0, s1, s3, t1, t2, t3, tUpd;
        auto huge = mVariables.mHuge;
        auto nGridX = mVariables.mGridX;
        auto nGridZ = mVariables.mGridZ;
        auto nCellX = mVariables.mCellX;
        auto nCellZ = mVariables.mCellZ;
        T t1, t2, t3, tUpd;
        T s0, s1, s3;
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        if (mVariables.mUniformGrid)
        {
            auto h = mVariables.mDx;
            for (int iz = nGridZ - 2; iz >= 0; --iz)
            {
                for (int ix = 1; ix < nGridX; ++ix)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    tUpd = cartesianFiniteDifference(huge, h,
                                                     s0, s1, s3,
                                                     t1, t2, t3);
                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        }
        else
        {
            auto dx = mVariables.mDx;
            auto dz = mVariables.mDz;
            auto dx_dz = mVariables.mDxDividedByDz;
            auto dz_dx = mVariables.mDzDividedByDx;
            auto cosTheta = mVariables.mCosTheta;
            auto sinTheta = mVariables.mSinTheta;
            for (int iz = nGridZ - 2; iz >= 0; --iz)
            {
                for (int ix=1; ix < nGridX; ++ix)
                {
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    tUpd = cartesianFiniteDifference(huge,
                                                     dx, dz,
                                                     dx_dz, dz_dx,
                                                     cosTheta, sinTheta,
                                                     s0, s1, s3,
                                                     t1, t2, t3);

                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        } // End check on uniform grid
    }
};

//----------------------------------------------------------------------------//
//                              Solves -x and -z                              //
//----------------------------------------------------------------------------//
template<class T>
class Solver2DSweep3 //: public Solver2DBaseClass<T>
{
public:
    Solver2DVariables<T> mVariables;
    /// Update travel time field with fast sweeping method
    void updateFSM(const T *slowness, T *travelTimes)
    {
        constexpr int sweep = 3;
        //T s0, s1, s3, t1, t2, t3, tUpd;
        auto huge = mVariables.mHuge;
        auto nGridX = mVariables.mGridX;
        auto nGridZ = mVariables.mGridZ;
        auto nCellX = mVariables.mCellX;
        auto nCellZ = mVariables.mCellZ;
        T t1, t2, t3, tUpd;
        T s0, s1, s3;
        int iCell0, iCell1, iCell3, it0, it1, it2, it3;
        if (mVariables.mUniformGrid)
        {
            auto h = mVariables.mDx;
            for (int iz = nGridZ - 2; iz >= 0; --iz)
            {
                for (int ix = nGridX - 2; ix >= 0; --ix)
                {
                    // Get surrounding slownesses
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];
                    // Get surrounding travel times
                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];
                    // Finite difference
                    tUpd = cartesianFiniteDifference(huge, h,
                                                     s0, s1, s3,
                                                     t1, t2, t3);
                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        }
        else
        {
            auto dx = mVariables.mDx;
            auto dz = mVariables.mDz;
            auto dx_dz = mVariables.mDxDividedByDz;
            auto dz_dx = mVariables.mDzDividedByDx;
            auto cosTheta = mVariables.mCosTheta;
            auto sinTheta = mVariables.mSinTheta;
            for (int iz = nGridZ - 2; iz >= 0; --iz)
            {
                for (int ix = nGridX - 2; ix >= 0; --ix)
                {
                    gridToSurroundingSlowness(sweep,
                                              ix, iz,
                                              nCellX, nCellZ,
                                              &iCell0, &iCell1, &iCell3);
                    s0 = slowness[iCell0];
                    s1 = slowness[iCell1];
                    s3 = slowness[iCell3];

                    gridToSurroundingTravelTimes(sweep, ix, iz, nGridX,
                                                 &it0, &it1, &it2, &it3);
                    t1 = travelTimes[it1];
                    t2 = travelTimes[it2];
                    t3 = travelTimes[it3];

                    tUpd = cartesianFiniteDifference(huge,
                                                     dx, dz,
                                                     dx_dz, dz_dx,
                                                     cosTheta, sinTheta,
                                                     s0, s1, s3,
                                                     t1, t2, t3);

                    travelTimes[it0] = std::min(travelTimes[it0], tUpd);
                }
            }
        } // End check on uniform grid
    }
};

/// Instantiate classes
//template class Solver2DBaseClass<double>;
//template class Solver2DBaseClass<float>;
template class Solver2DSweep0<double>;
template class Solver2DSweep0<float>;
template class Solver2DSweep1<double>;
template class Solver2DSweep1<float>;
template class Solver2DSweep2<double>;
template class Solver2DSweep2<float>;
template class Solver2DSweep3<double>;
template class Solver2DSweep3<float>;
*/

template class Solver2DSweep<double, EikonalXX::SweepNumber2D::SWEEP1>;
template class Solver2DSweep<double, EikonalXX::SweepNumber2D::SWEEP2>;
template class Solver2DSweep<double, EikonalXX::SweepNumber2D::SWEEP3>;
template class Solver2DSweep<double, EikonalXX::SweepNumber2D::SWEEP4>;
template class Solver2DSweep<float, EikonalXX::SweepNumber2D::SWEEP1>;
template class Solver2DSweep<float, EikonalXX::SweepNumber2D::SWEEP2>;
template class Solver2DSweep<float, EikonalXX::SweepNumber2D::SWEEP3>;
template class Solver2DSweep<float, EikonalXX::SweepNumber2D::SWEEP4>;
}
#endif
