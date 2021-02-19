#include <iostream>
#include <vector>
#include <algorithm>
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
    std::vector<int> mLevelPointer;
    /// Map from a node to the corresponding level.  This has dimension
    /// mGrid. 
    std::vector<int> mNodeToLevel;
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
    pImpl->mLevelPointer.clear();
    pImpl->mNodeToLevel.clear();
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
    // Compute level pointer
    std::vector<int> work(pImpl->mGrid, 0);
    pImpl->mLevelPointer.resize(pImpl->mNumberOfLevels + 1, -1);
    auto levelPtr = pImpl->mLevelPointer.data();
    levelPtr[0] = 0;
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
        for (int iz = k1; iz <= k2; ++iz)
        {
            auto j1 = std::max(0, level - iz - (nx - 1));
            j1 = std::min(j1, ny - 1); // Don't include boundary layer
            auto j2 = std::min(ny - 1, level - iz);
            for (int iy = j1; iy <= j2; ++iy)
            {
                auto i1 = std::max(0, level - iz - iy);
                i1 = std::min(i1, nx - 1); // Don't include boundary layer
                auto i2 = std::min(nx - 1, level - iz - iy);
                for (int ix = i1; ix <= i2; ++ix)
                {
std::cout << level << " " << ix << " " << iy << " " << iz << " " << gridToIndex(nx, ny, ix, iy, iz) << std::endl;
                    work[nNodesInLevel] = gridToIndex(nx, ny, ix, iy, iz);
                    nNodesInLevel = nNodesInLevel + 1; 
                }
            } 
        }
std::cout << std::endl;
        pImpl->mMaxLevelSize = std::max(pImpl->mMaxLevelSize, nNodesInLevel);
        levelPtr[level + 1] = levelPtr[level] + nNodesInLevel;
        std::sort(work.data(), work.data() + nNodesInLevel);
        // Convert grid to node
        for (int i = 0; i <nNodesInLevel; ++i)
        {
            int ix, iy, iz;
            indexToGrid(work[i], nx, ny, &ix, &iy, &iz);
        }
    }
    pImpl->mInitialized = true;
}

/// Class initialized?
template<EikonalXX::SweepNumber3D E>
bool Graph3D<E>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

template<EikonalXX::SweepNumber3D E>
int Graph3D<E>::getNumberOfLevels() const
{
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    return pImpl->mNumberOfLevels;
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
