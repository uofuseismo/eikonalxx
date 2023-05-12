#ifndef EIKONALXX_PRIVATE_SOLVER3D_HPP
#define EIKONALXX_PRIVATE_SOLVER3D_HPP
#include <iomanip>
#include <limits>
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/graph3d.hpp"
#include "solverUtilities3d.hpp"


#define HUGE 1.e10

namespace
{

using namespace EikonalXX;

/// Solver class
template<class T, EikonalXX::SweepNumber3D E>
class Solver3DSweep
{
public:
    /// Initialize the variables in this class 
    void initialize(const EikonalXX::SolverOptions &options,
                    const EikonalXX::Geometry3D &geometry)
    {
        // Solver options
        mAlgorithm = options.getAlgorithm();
        mFactoredEikonalSolverRadius
            = options.getFactoredEikonalEquationSolverRadius(); 
        mConvergenceTolerance = options.getConvergenceTolerance();
        // Uniform grid?
        auto dx = geometry.getGridSpacingInX();
        auto dy = geometry.getGridSpacingInY();
        auto dz = geometry.getGridSpacingInZ();
        mUniformGrid = false;
        if (std::abs(dx - dy) < std::numeric_limits<T>::epsilon()*100 &&
            std::abs(dx - dz) < std::numeric_limits<T>::epsilon()*100)
        {
            mUniformGrid = true;
        }
        // Define the grid spacing.
        mDx = static_cast<T> (dx);
        mDy = static_cast<T> (dy);
        mDz = static_cast<T> (dz);
        // Define the grid size.
        mGridX = geometry.getNumberOfGridPointsInX();
        mGridY = geometry.getNumberOfGridPointsInY();
        mGridZ = geometry.getNumberOfGridPointsInZ();
        mCellX = geometry.getNumberOfCellsInX();
        mCellY = geometry.getNumberOfCellsInY();
        mCellZ = geometry.getNumberOfCellsInZ();
        // Initialize graph and allocate space for level set method
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LevelSetMethod)
        {
            // Create graph
            mGraph.initialize(mGridX, mGridY, mGridZ);
            auto nLevels = mGraph.getNumberOfLevels();
            // Set space for sweep's slowness
            mSweepSlowness.resize(nLevels);
            for (int level = 0; level < nLevels; ++level)
            {
                auto nNodesInLevel = mGraph.getNumberOfNodesInLevel(level);
                mSweepSlowness[level].allocate(nNodesInLevel);
            }
        }

    }
    /// Set the velocity model
    void setVelocityModel(const EikonalXX::Model3D<T> &velocityModel)
    {
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LevelSetMethod)
        {
#ifndef NDEBUG
            for (int level = 0; level < mLevels; ++level)
            {
                mSweepSlowness[level].zero();
            }
#endif
            // Fill slownesses in each sweep
            slownessToSweepSlowness<T, E>(mCellX, mCellY, mCellZ,
                                          mGraph,
                                          velocityModel.getSlownessPointer(),
                                          mSweepSlowness.data());
#ifndef NDEBUG
            for (int level = 0; level < mLevels; ++level)
            {
                assert(mSweepSlowness[level].getMinimumValue() > 0); 
            }
#endif
        }
    } 
    /// Clears the class / releases memory
    void clear() noexcept
    {
        mSweepSlowness.clear();
        mGraph.clear();
        mDx = 0;
        mDy = 0;
        mDz = 0;
        mMaxTravelTimeChange = HUGE;
        mSourceSlowness = 0;
        mConvergenceTolerance = 0;
        mSourceOffsetX = 0;
        mSourceOffsetY = 0;
        mSourceOffsetZ = 0;
        mGridX = 0;
        mGridY = 0;
        mGridZ = 0;
        mCellX = 0;
        mCellY = 0;
        mCellZ = 0;
        mSourceIndexX = 0;
        mSourceIndexY = 0;
        mSourceIndexZ = 0;
        mLevels = 0;
        mFactoredEikonalSolverRadius = 0;
        mAlgorithm = EikonalXX::SolverAlgorithm::LevelSetMethod;
        mUniformGrid = true;
    }
    /// Sets the source information
    /// @param[in] iSrcX  The source cell in x.
    /// @param[in] iSrcY  The source cell in y.
    /// @param[in] iSrcZ  The source cell in z.
    /// @param[in] xSourceOffset   The x source position in meters from the
    ///                            x origin.
    /// @param[in] ySourceOffset   The y source position in meters from the
    ///                            y origin.
    /// @param[in] zSourceOffset   The z source position in meters from the
    ///                            z origin.
    /// @param[in] sourceSlowness  The slowness in s/m in the source cell. 
    void setSourceInformation(const int iSrcX, const int iSrcY, const int iSrcZ,
                              const T xSourceOffset,
                              const T ySourceOffset,
                              const T zSourceOffset,
                              const T sourceSlowness)
    {   
        mSourceIndexX = iSrcX;
        mSourceIndexY = iSrcY;
        mSourceIndexZ = iSrcZ;
        mSourceOffsetX = xSourceOffset;
        mSourceOffsetY = ySourceOffset;
        mSourceOffsetZ = zSourceOffset;
        mSourceSlowness = sourceSlowness;
    }
    /// Performs the fast sweeping method update.
    /// @param[in] slowness         The cell-based slowness model in s/m.
    /// @param[in,out] travelTimes  The traveltimes at all grid points in s.
    /// @param[in] initialize       True indicates that this is an
    ///                             initialization sweep.
    void updateFSM(const T *__restrict__ slowness, T *__restrict__ travelTimes,
                   const bool initialize)
    {
std::cout.precision(16);
        // Get the loop limits
        int ix0, ix1, ixDir, iy0, iy1, iyDir, iz0, iz1, izDir;
        if (initialize)
        {
            getLoopLimits<E>(mGridX, mGridY, mGridZ,
                             mSourceIndexX, mSourceIndexY, mSourceIndexZ,
                             &ix0, &iy0, &iz0,
                             &ix1, &iy1, &iz1, 
                             &ixDir, &iyDir, &izDir);
        }
        else
        {
            getLoopLimits<E>(mGridX, mGridY, mGridZ,
                             &ix0, &iy0, &iz0,
                             &ix1, &iy1, &iz1, 
                             &ixDir, &iyDir, &izDir);
        }
        // Get finite difference sweep signs
        int ixShift, iyShift, izShift, signX, signY, signZ;
        getSweepFiniteDifferenceSigns<E>(&ixShift, &iyShift, &izShift,
                                         &signX, &signY, &signZ);
        // Some local variables
        constexpr T huge = HUGE;
        T t0, t1, t2, t3, t4, t5, t6, t7, tUpd;
        T s0, s1, s2, s3, s4, s5, s7; 
        int iCell0, iCell1, iCell2, iCell3, iCell4, iCell5, iCell7;
        int it0, it1, it2, it3, it4, it5, it6, it7;
        T diffMax = HUGE;
        // Uniform finite difference stencil 
        if (mUniformGrid)
        {
            auto h = mDx; // dx = dy = dz
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int iy = iy0; iy  != iy1; iy = iy + iyDir)
                {
                    for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                    {
                        // Get surrounding slownesses
                        gridToSurroundingSlownessIndices<E>(
                            ix, iy, iz, 
                            mCellX, mCellY, mCellZ,
                            &iCell0, &iCell1, &iCell2, &iCell3,
                            &iCell4, &iCell5, &iCell7);
                        s0 = slowness[iCell0];
                        s1 = slowness[iCell1];
                        s2 = slowness[iCell2]; 
                        s3 = slowness[iCell3];
                        s4 = slowness[iCell4];
                        s5 = slowness[iCell5];
                        s7 = slowness[iCell7];
//std::cout << s0 << " " << s1 << " " << s2 << " " << s3 << " " << s4 << " " << s5 << " " << s7 << std::endl;
                        // Get surrounding travel times
                        gridToSurroundingTravelTimeIndices<E>(
                            ix, iy, iz,
                            mGridX, mGridY, mGridZ,
                            &it0, &it1, &it2, &it3,
                            &it4, &it5, &it6, &it7);
                        t0 = travelTimes[it0];
                        t1 = travelTimes[it1];
                        t2 = travelTimes[it2];
                        t3 = travelTimes[it3];
                        t4 = travelTimes[it4];
                        t5 = travelTimes[it5];
                        t6 = travelTimes[it6];
                        t7 = travelTimes[it7];
//std::cout << t1 << " " << t2 << " " << t3 << " " << t3 << " " << t4 << " " << t5 << " " << t6 << " " << t7 << std::endl;
                        // Finite difference
                        tUpd = ::finiteDifference(mFactoredEikonalSolverRadius,
                                 huge,
                                 h, 
                                 mSourceSlowness,
                                 ix, iy, iz,
                                 signX, signY, signZ,
                                 ixShift, iyShift, izShift,
                                 mSourceIndexX, mSourceIndexY, mSourceIndexZ,
                                 mSourceOffsetX, mSourceOffsetY, mSourceOffsetZ,
                                 s0, s1, s2, s3,
                                 s4, s5, s7,
                                 t1, t2, t3,
                                 t4, t5, t6, t7);
                        // Update node if new travel time is smaller
                        if (tUpd < t0)
                        {
                            travelTimes[it0] = tUpd;
                            diffMax = sycl::fmax(t0 - tUpd, diffMax);
                        }
                    }
                }
            }
        }
        else
        {
            T dxInv = 1/mDx;
            T dyInv = 1/mDy;
            T dzInv = 1/mDz;
            T dx2Inv = 1/(mDx*mDx);
            T dy2Inv = 1/(mDy*mDy);
            T dz2Inv = 1/(mDz*mDz);
            T dx_dz = mDx/mDz;
            T dz_dx = mDz/mDx;
            T dy_dz = mDy/mDz;
            T dz_dy = mDz/mDy;
            T dx_dy = mDx/mDy;
            T dy_dx = mDy/mDx;
            T dx2_p_dz2_inv = 1/(mDx*mDx + mDz*mDz);
            T dy2_p_dz2_inv = 1/(mDy*mDy + mDz*mDz);
            T dx2_p_dy2_inv = 1/(mDx*mDx + mDy*mDy);
            for (int iz = iz0; iz != iz1; iz = iz + izDir)
            {
                for (int iy = iy0; iy  != iy1; iy = iy + iyDir)
                {
                    for (int ix = ix0; ix != ix1; ix = ix + ixDir)
                    {
                        gridToSurroundingSlownessIndices<E>(
                            ix, iy, iz, 
                            mCellX, mCellY, mCellZ,
                            &iCell0, &iCell1, &iCell2, &iCell3,
                            &iCell4, &iCell5, &iCell7);
                        s0 = slowness[iCell0];
                        s1 = slowness[iCell1];
                        s2 = slowness[iCell2];
                        s3 = slowness[iCell3];
                        s4 = slowness[iCell4];
                        s5 = slowness[iCell5];
                        s7 = slowness[iCell7];

                        gridToSurroundingTravelTimeIndices<E>(
                            ix, iy, iz,
                            mGridX, mGridY, mGridZ,
                            &it0, &it1, &it2, &it3,
                            &it4, &it5, &it6, &it7);
                        t0 = travelTimes[it0];
                        t1 = travelTimes[it1];
                        t2 = travelTimes[it2];
                        t3 = travelTimes[it3];
                        t4 = travelTimes[it4];
                        t5 = travelTimes[it5];
                        t6 = travelTimes[it6];
                        t7 = travelTimes[it7];

                        tUpd = ::finiteDifference(
                                 mFactoredEikonalSolverRadius,
                                 huge,
                                 mDx, mDy, mDz,
                                 dxInv, dyInv, dzInv,
                                 dx2Inv, dy2Inv, dz2Inv,
                                 dx_dz, dz_dx,
                                 dy_dz, dz_dy,
                                 dy_dx, dx_dy,
                                 dx2_p_dz2_inv, dy2_p_dz2_inv, dx2_p_dy2_inv,
                                 mSourceSlowness,
                                 ix, iy, iz,
                                 signX, signY, signZ,
                                 ixShift, iyShift, izShift, 
                                 mSourceIndexX, mSourceIndexY, mSourceIndexZ,
                                 mSourceOffsetX, mSourceOffsetY, mSourceOffsetZ,
                                 s0, s1, s2, s3,
                                 s4, s5, s7,
                                 t1, t2, t3,
                                 t4, t5, t6, t7);

                        if (tUpd < t0)
                        {
                            travelTimes[it0] = tUpd;
                            diffMax = sycl::fmax(t0 - tUpd, diffMax);
                        }
                    }
                }
            }
        }
        mMaxTravelTimeChange = diffMax;
    }
//private:
    /// Holds the slowness in the sweep.
    std::vector<SweepSlowness3D<T>> mSweepSlowness;
    /// Defines the level set method elimination graph.
    Graph3D<E> mGraph;
    /// Grid spacing in x, y, and z.
    T mDx = 0;
    T mDy = 0;
    T mDz = 0;
    /// The source offset in x, y, and z.
    T mSourceOffsetX = 0;
    T mSourceOffsetY = 0;
    T mSourceOffsetZ = 0;
    /// The slowness at the source in s/m.
    T mSourceSlowness = 0;
    /// The largest change in the travel time in s for this sweep.
    T mMaxTravelTimeChange = HUGE;
    /// The convergence tolerance in seconds.
    T mConvergenceTolerance = 0;
    /// Number of grid points in x, y, and z.
    int mGridX = 0;
    int mGridY = 0;
    int mGridZ = 0;
    /// Number of cells in x, y, and z.
    int mCellX = 0;
    int mCellY = 0;
    int mCellZ = 0;
    /// The source grid indices.
    int mSourceIndexX = 0;
    int mSourceIndexY = 0;
    int mSourceIndexZ = 0;
    /// Number of levels.
    int mLevels = 0;
    /// Spherical to cartesian transition.
    int mFactoredEikonalSolverRadius = 0;
    /// Algorithm type.
    EikonalXX::SolverAlgorithm mAlgorithm
        = EikonalXX::SolverAlgorithm::LevelSetMethod;
    /// Is the grid spacing the same in x, y, and z?
    bool mUniformGrid = true;
};

///--------------------------------------------------------------------------///
///                               Template Instantiation                     ///
///--------------------------------------------------------------------------///

template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep1>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep2>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep3>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep4>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep5>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep6>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep7>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::Sweep8>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep1>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep2>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep3>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep4>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep5>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep6>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep7>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::Sweep8>;


}
#endif
