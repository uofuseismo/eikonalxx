#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#ifdef USE_BOOST
//#include <boost/graph/properties.hpp>
//#include <boost/graph/topological_sort.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/edge_coloring.hpp>
#endif
#include "eikonalxx/graph3d.hpp"
#include "private/solverUtilities3d.hpp"

using namespace EikonalXX;

template<EikonalXX::SweepNumber3D E>
class Graph3D<E>::Graph3DImpl
{
public:
    /// Map from the level to the start index of the nodes in that level.
    /// This has dimension [mNumberOfLevels + 1].
    std::vector<int> mLevelStartPointer;
    /// Map from a global grid idnex to the corresponding level.  This has
    /// dimension [mGrid]. 
    std::vector<int> mIndexToLevel;
    /// Converts the node in a level to the global grid index.
    /// This is an array whose dimension is [mGrid] but is accessed with
    /// mLevelStartPointer.
    std::vector<int> mNodeInLevelToIndex;
    /// Converts the node in a level to the corresponding grid point in x.
    /// This is an array whose dimension is [mGrid] but is accessed with
    /// mLevelStartPointer.
    std::vector<int> mNodeInLevelToXGridPoint;
    /// Converts the node in a level to the corresponding grid point in y.
    /// This is an array whose dimension is [mGrid] but is accessed with
    /// mLevelStartPointer.
    std::vector<int> mNodeInLevelToYGridPoint;
    /// Converts the node in a level to the corresponding grid point in z.
    /// This is an array whose dimension is [mGrid] but is accessed with
    /// mLevelStartPointer.
    std::vector<int> mNodeInLevelToZGridPoint;
    /// Number of grid points in (x,y,z)
    int mNx = 0;
    int mNy = 0;
    int mNz = 0;
    /// Number of grid points (= mNx*mNy*mNz)
    int mGrid = 0;
    /// Number of levels
    int mNumberOfLevels = 0;
    /// Maximum size of a level
    int mMaxLevelSize = 0;
    /// Class initialized?
    bool mInitialized = false;
};

/// C'tor
template<EikonalXX::SweepNumber3D E>
Graph3D<E>::Graph3D() :
    pImpl(std::make_unique<Graph3DImpl> ())
{
}

/// Copy c'tor
template<EikonalXX::SweepNumber3D E>
Graph3D<E>::Graph3D(const Graph3D &graph)
{
    *this = graph;
}

/// Move c'tor
template<EikonalXX::SweepNumber3D E>
Graph3D<E>::Graph3D(Graph3D &&graph) noexcept
{
    *this = std::move(graph);
}

/// Releases memory and resets
template<EikonalXX::SweepNumber3D E>
void Graph3D<E>::clear() noexcept
{
    pImpl->mLevelStartPointer.clear();
    pImpl->mNodeInLevelToIndex.clear();
    pImpl->mNodeInLevelToXGridPoint.clear();
    pImpl->mNodeInLevelToYGridPoint.clear();
    pImpl->mNodeInLevelToZGridPoint.clear();
    pImpl->mIndexToLevel.clear();
    pImpl->mNx = 0;
    pImpl->mNy = 0;
    pImpl->mNz = 0;
    pImpl->mGrid = 0;
    pImpl->mNumberOfLevels = 0;
    pImpl->mMaxLevelSize = 0;
    pImpl->mInitialized = false;
}

/// Destructor
template<EikonalXX::SweepNumber3D E>
Graph3D<E>::~Graph3D() = default;

/// Copy assignment
template<EikonalXX::SweepNumber3D E>
Graph3D<E>& Graph3D<E>::operator=(const Graph3D<E> &graph)
{
    if (&graph == this){return *this;}
    pImpl = std::make_unique<Graph3DImpl> (*graph.pImpl);
    return *this;
}

/// Move assignment
template<EikonalXX::SweepNumber3D E>
Graph3D<E>& Graph3D<E>::operator=(Graph3D<E> &&graph) noexcept
{
    if (&graph == this){return *this;}
    pImpl = std::move(graph.pImpl);
    return *this;
}

/// Initialize
template<EikonalXX::SweepNumber3D E>
void Graph3D<E>::initialize(const int nx, const int ny, const int nz)
{
    clear();
    if (nx < 3){throw std::invalid_argument("nx must be at least 3");}
    if (ny < 3){throw std::invalid_argument("ny must be at least 3");}
    if (nz < 3){throw std::invalid_argument("nz must be at least 3");}
    pImpl->mNumberOfLevels = EikonalXX::computeNumberOfLevels(nx, ny, nz);
    pImpl->mNx = nx;
    pImpl->mNy = ny;
    pImpl->mNz = nz;
    pImpl->mGrid = pImpl->mNx*pImpl->mNy*pImpl->mNz;
    // Allocate space
    pImpl->mIndexToLevel.resize(pImpl->mGrid, -1);
    pImpl->mNodeInLevelToIndex.resize(pImpl->mGrid, -1);
    pImpl->mNodeInLevelToXGridPoint.resize(pImpl->mGrid, -1);
    pImpl->mNodeInLevelToYGridPoint.resize(pImpl->mGrid, -1);
    pImpl->mNodeInLevelToZGridPoint.resize(pImpl->mGrid, -1);
    // Get pointers
    pImpl->mLevelStartPointer.resize(pImpl->mNumberOfLevels + 1, -1);
    auto indexToLevelPtr = pImpl->mIndexToLevel.data();
    auto nodeInLevelToIndexPtr = pImpl->mNodeInLevelToIndex.data();
    auto nodeInLevelToXGridPointPtr = pImpl->mNodeInLevelToXGridPoint.data();
    auto nodeInLevelToYGridPointPtr = pImpl->mNodeInLevelToYGridPoint.data();
    auto nodeInLevelToZGridPointPtr = pImpl->mNodeInLevelToZGridPoint.data();
    auto levelStartPtr = pImpl->mLevelStartPointer.data();
    // Build up pointers
    levelStartPtr[0] = 0;
    for (int level = 0; level < pImpl->mNumberOfLevels; ++level)
    {
        // Solve level = constant = ix + iy + iz.
        // The idea to find the minimum z value is to solve for the maximum
        // allowable nx and ny values.  Likewise, the end value should be
        // the minimum allowable x and z values which are 0.  Hence,
        // for k2 we simplify z = level - x - y = level - 0 - 0 = level.
        // The idea is then repeated for y where we choose the specific
        // value of iz.
        auto k1 = std::max(0, level - (nx - 1) - (ny - 1));
        k1 = std::min(k1, nz - 1);
        auto k2 = std::min(nz - 1, level);
        int nNodesInLevel = 0;
        int offset = levelStartPtr[level];
        for (int iz = k1; iz <= k2; ++iz)
        {
            auto j1 = std::max(0, level - iz - (nx - 1));
            j1 = std::min(j1, ny - 1); // Ensure y start is in bounds
            auto j2 = std::min(ny - 1, level - iz);
            for (int iy = j1; iy <= j2; ++iy)
            {
                auto i1 = std::max(0, level - iz - iy);
                i1 = std::min(i1, nx - 1); // Ensure x start is in bounds 
                auto i2 = std::min(nx - 1, level - iz - iy);
                for (int ix = i1; ix <= i2; ++ix)
                {
                    int jx, jy, jz;
                    permuteGrid<E>(ix, iy, iz, nx, ny, nz, &jx, &jy, &jz);
                    nodeInLevelToIndexPtr[offset + nNodesInLevel]
                        = gridToIndex(nx, ny, jx, jy, jz);
                    indexToLevelPtr[offset + nNodesInLevel] = level;
                    nNodesInLevel = nNodesInLevel + 1; 
                }
            } 
        }
        // Update pointers
        pImpl->mMaxLevelSize = std::max(pImpl->mMaxLevelSize, nNodesInLevel);
        levelStartPtr[level + 1] = levelStartPtr[level] + nNodesInLevel;
        // Sort in increasing grid order to make travel time table accesses in
        // solver phase slightly less random and promote a little cache
        // coherency.
        std::sort(nodeInLevelToIndexPtr + offset,
                  nodeInLevelToIndexPtr + offset + nNodesInLevel);
        // Convert grid indices in level to (ix, iy, iz) grid point
        #pragma omp simd
        for (int i = 0; i < nNodesInLevel; ++i)
        {
            int ix, iy, iz;
            indexToGrid(nodeInLevelToIndexPtr[offset + i],
                        nx, ny, &ix, &iy, &iz);
            //permuteGrid<E>(ix, iy, iz, nx, ny, nz, &jx, &jy, &jz);
//std::cout << static_cast<int> (E) << " " << level << " " << ix << " " << iy << " " << iz << " " << nodeInLevelToIndexPtr[offset + i] << std::endl;
            nodeInLevelToXGridPointPtr[offset + i] = ix;
            nodeInLevelToYGridPointPtr[offset + i] = iy;
            nodeInLevelToZGridPointPtr[offset + i] = iz;
        }
    }
#ifndef NDEBUG
    auto mm = std::minmax_element(pImpl->mIndexToLevel.begin(),
                                  pImpl->mIndexToLevel.end());
    assert(*mm.first == 0);
    assert(*mm.second == pImpl->mNumberOfLevels - 1);
    mm = std::minmax_element(pImpl->mNodeInLevelToXGridPoint.begin(),
                             pImpl->mNodeInLevelToXGridPoint.end());
    assert(*mm.first == 0);
    assert(*mm.second == pImpl->mNx - 1);
    mm = std::minmax_element(pImpl->mNodeInLevelToYGridPoint.begin(),
                             pImpl->mNodeInLevelToYGridPoint.end());
    assert(*mm.first == 0);
    assert(*mm.second == pImpl->mNy - 1);
    mm = std::minmax_element(pImpl->mNodeInLevelToZGridPoint.begin(),
                             pImpl->mNodeInLevelToZGridPoint.end());
    assert(*mm.first == 0);
    assert(*mm.second == pImpl->mNz - 1);
#endif
    pImpl->mInitialized = true;
}

/// Class initialized?
template<EikonalXX::SweepNumber3D E>
bool Graph3D<E>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Number of grid points
template<EikonalXX::SweepNumber3D E>
int Graph3D<E>::getNumberOfGridPoints() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mGrid;
}

/// Number of levels
template<EikonalXX::SweepNumber3D E>
int Graph3D<E>::getNumberOfLevels() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mNumberOfLevels;
}

/// Max level size
template<EikonalXX::SweepNumber3D E>
int Graph3D<E>::getMaximumLevelSize() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mMaxLevelSize;
}

/// Map from grid index to level
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getIndexToLevelPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mIndexToLevel.data();
}

/// Map from node in level to global index
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getNodeInLevelToIndexPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mNodeInLevelToIndex.data();
}

/// Map from level to start index in nodeInLevel arrays
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getLevelStartPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mLevelStartPointer.data();
}

/// Map from the node in level to x grid point
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getNodeInLevelToXGridPointPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mNodeInLevelToXGridPoint.data();
}

/// Map from the node in level to y grid point
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getNodeInLevelToYGridPointPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mNodeInLevelToYGridPoint.data();
}

/// Map from the node in level to z grid point
template<EikonalXX::SweepNumber3D E>
const int *Graph3D<E>::getNodeInLevelToZGridPointPointer() const
{
    if (!isInitialized()){throw std::runtime_error("Class is not initialized");}
    return pImpl->mNodeInLevelToZGridPoint.data();
}

/// @brief From the number of grid points this returns the number of levels.
/// @param[in] nx  The number of grid points in x.
/// @param[in] ny  The number of grid points in y.
/// @param[in] nz  The number of grid points in z.
/// @result The number of levels in the 3D solver.
int EikonalXX::computeNumberOfLevels(const int nx,
                                     const int ny,
                                     const int nz) noexcept
{
    return nx + ny + nz - 2;
}
///--------------------------------------------------------------------------///
///                           Template Instantiation                         ///
///--------------------------------------------------------------------------///
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP1>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP2>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP3>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP4>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP5>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP6>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP7>;
template class EikonalXX::Graph3D<EikonalXX::SweepNumber3D::SWEEP8>;
