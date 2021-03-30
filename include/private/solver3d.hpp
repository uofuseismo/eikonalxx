#ifndef EIKONALXX_PRIVATE_SOLVER3D_HPP
#define EIKONALXX_PRIVATE_SOLVER3D_HPP
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
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD)
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
        if (mAlgorithm == EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD)
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
        mGridX = 0;
        mGridY = 0;
        mGridZ = 0;
        mCellX = 0;
        mCellY = 0;
        mCellZ = 0;
        mSphericalSolverRadius = 0;
        mAlgorithm = EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD;
        mUniformGrid = true;
    }
    /// Performs the fast sweeping method on a grid
    void updateFSM(const T *slowness, T *travelTimes)
    {   
        int ix0, ix1, ixDir, iy0, iy1, iyDir, iz0, iz1, izDir;
        getLoopLimits<E>(mGridX, mGridY, mGridZ,
                         &ix0, &iy0, &iz0,
                         &ix1, &iy1, &iz1, 
                         &ixDir, &iyDir, &izDir);
/*
        int ixShift, izShift, signX, signZ;
        getSweepFiniteDifferenceSigns(E,
                                      &ixShift, &izShift,
                                      &signX, &signZ);
*/
        T t1, t2, t3, t4, t5, t6, t7, tUpd;
        T s0, s1, s2, s3, s4, s5, s7; 
        int iCell0, iCell1, iCell2, iCell3, iCell4, iCell5, iCell7;
        int it0, it1, it2, it3, it4, it5, it6, it7;
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
                        // Get surrounding travel times
                        gridToSurroundingTravelTimeIndices<E>(
                            ix, iy, iz,
                            mGridX, mGridY, mGridZ,
                            &it0, &it1, &it2, &it3,
                            &it4, &it5, &it6, &it7);
                    }
                }
            }
        }
        else
        {
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
                    }
                }
            }
        }
    }
//private:
    /// Holds the slowness in the sweep.
    std::vector<SweepSlowness3D<T>> mSweepSlowness;
    /// Defines the level set method elimination graph.
    Graph3D<E> mGraph;
    /// Grid spacing in x, y, and z.
    float mDx = 0;
    float mDy = 0;
    float mDz = 0;
    /// Number of grid points in x, y, and z.
    int mGridX = 0;
    int mGridY = 0;
    int mGridZ = 0;
    /// Number of cells in x, y, and z.
    int mCellX = 0;
    int mCellY = 0;
    int mCellZ = 0;
    /// Number of levels.
    int mLevels = 0;
    /// Spherical to cartesian transition.
    int mSphericalSolverRadius = 0;
    /// Algorithm type.
    EikonalXX::SolverAlgorithm mAlgorithm
        = EikonalXX::SolverAlgorithm::LEVEL_SET_METHOD;
    /// Is the grid spacing the same in x, y, and z?
    bool mUniformGrid = true;
};

///--------------------------------------------------------------------------///
///                               Template Instantiation                     ///
///--------------------------------------------------------------------------///

template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP1>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP2>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP3>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP4>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP5>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP6>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP7>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP8>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP1>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP2>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP3>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP4>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP5>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP6>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP7>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP8>;


}
#endif
