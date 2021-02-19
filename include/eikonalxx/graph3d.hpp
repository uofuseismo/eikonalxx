#ifndef EIKONALXX_GRAPH3D_HPP
#define EIKONALXX_GRAPH3D_HPP
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
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

    /// @name Destructors
    /// @{
    /// @brief Destructor.
    ~Graph3D();
    /// @brief Releases the memory on the class.
    void clear() noexcept;
    /// @}

    /// @result The maximum number of nodes in a level.
    /// @throws std::invalid_argument if \c isInitialized() is false.
    int getMaximumLevelSize() const;
    /// @result The number of levels.
    /// @throws std::invalid_argument if \c isInitialized() is false.
    int getNumberOfLevels() const;
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
