#ifndef EIKONALXX_GRAPH_3D_HPP
#define EIKONALXX_GRAPH_3D_HPP
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
/// @class Graph3D "graph3d.hpp" "eikonalxx/graph3d.hpp"
/// @brief Defines the computational graph for the 3D level-set method.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<EikonalXX::SweepNumber3D E>
class Graph3D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Graph3D();
    /// @brief Copy constructor.
    /// @param[in] graph  The graph class from which to initialize this class.
    Graph3D(const Graph3D &graph);
    /// @brief Move constructor.
    /// @param[in,out] graph  The graph class from which to initialize this
    ///                       class.  On exit, graph's behavior is undefined.
    Graph3D(Graph3D &&graph) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] graph   The graph to copy to this.
    /// @result A deep copy of the input graph.
    Graph3D& operator=(const Graph3D &graph);
    /// @brief Move assignment operator.
    /// @param[in,out] graph  The graph whose memory will be moved to this.
    ///                       On exit, graph's behavior is undefined. 
    /// @result The memory from graph moved to this.
    Graph3D& operator=(Graph3D &&graph) noexcept;
    /// @}

    /// @name Initialization
    /// @{

    /// @brief Initializes the graph.
    /// @param[in] nx   The number of grid points in x.
    /// @param[in] ny   The number of grid points in y.
    /// @param[in] nz   The number of grid points in z.
    /// @throws std::invalid_argument if nx, ny, or nz is not at least 3.
    void initialize(int nx, int ny, int nz);
    /// @result True indicates that the class is initialized.
    bool isInitialized() const noexcept;
    /// @}

    /// @result The number of grid points in the model.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfGridPoints() const;
    /// @result The maximum number of nodes in a level.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getMaximumLevelSize() const;
    /// @result The number of levels.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfLevels() const;
    /// @param[in] level  The level for which the number of nodes is requested.
    /// @result The number of nodes in the level'th level.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if level is not in the range of
    ///         [0, \c getNumberOfLevels() - 1].
    [[nodiscard]] int getNumberOfNodesInLevel(int level) const;
    /// @result A map from a global grid index to the level.  This is an array
    ///         whose dimension is [\c getNumberOfGridPoints()].
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getIndexToLevelPointer() const;
    /// @result The level'th index maps to the start index of the
    ///         nodes in the level'th level.  The number of nodes in
    ///         the level'th level are given by
    ///         levelPtr[level+1] - levelPtr[level].  This is an array
    ///         whose dimension is [\c getNumberOfLevels() + 1].
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getLevelStartPointer() const;
    /// @result This is a map from the node in a level to the global grid
    ///         index value.  This is an array whose dimension is 
    ///         [ \c getNumberOfGridPoints() ].
    ///         @sa \c getLevelStartPointer().
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getNodeInLevelToIndexPointer() const;
    /// @result This is a map from the node in a level to the x grid
    ///         point.  This is an array whose dimension is
    ///         [ \c getNumberOfGridPoints() ]. 
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getNodeInLevelToXGridPointPointer() const;
    /// @result This is a map from the node in a level to the y grid
    ///         point.  This is an array whose dimension is
    ///         [ \c getNumberOfGridPoints() ]. 
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getNodeInLevelToYGridPointPointer() const;
    /// @result This is a map from the node in a level to the z grid
    ///         point.  This is an array whose dimension is
    ///         [ \c getNumberOfGridPoints() ]. 
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] const int *getNodeInLevelToZGridPointPointer() const;

    /// @name Destructors
    /// @{

    /// @brief Releases the memory on the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Graph3D();
    /// @}
private:
    class Graph3DImpl;
    std::unique_ptr<Graph3DImpl> pImpl;
};
/// @brief From the number of grid points this returns the number of levels.
/// @param[in] nx  The number of grid points in x.
/// @param[in] ny  The number of grid points in y.
/// @param[in] nz  The number of grid points in z.
/// @result The number of levels in the 3D solver.
int computeNumberOfLevels(const int nx, const int ny, const int nz) noexcept;
}
#endif
