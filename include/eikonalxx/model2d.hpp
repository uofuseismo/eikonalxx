#ifndef EIKONALXX_MODEL2D_HPP
#define EIKONALXX_MODEL2D_HPP
#include <memory>
#include <vector>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
// Forward declarations
class Geometry2D;
/*!
 * @class Model2D "model2d.hpp" "eikonalxx/model2d.hpp"
 * @brief Defines a 2D velocity model.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
template<class T>
class Model2D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Default constructor.
    Model2D();
    /// @brief Copy constructor.
    /// @param[in] model  The model class from which to initialize this class.
    Model2D(const Model2D &model);
    /// @brief Move constructor.
    /// @param[in,out] model  The model class from which to initialize this
    ///                       class.  On exit, model's behavior is undefined.
    Model2D(Model2D &&model) noexcept;
    /// @} 

    /// @name Operators
    /// @{
    /// @brief Copy assignment.
    /// @param[in] model  The model to copy to this.
    /// @result A deep copy of model to this.
    Model2D& operator=(const Model2D &model);
    /// @brief Move assignment.
    /// @param[in,out] model  The memory on model to move to this.
    /// @result The memory from model moved to this.
    Model2D& operator=(Model2D &&model) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Destructor.
    ~Model2D();
    /// @brief Resets the class and releases all memory.
    void clear() noexcept;
    /// @}

    /// @brief Initializes the class.
    /// @param[in] geometry   The geometry class.  This must have the number of
    ///                       grid points in x and z as well as the grid
    ///                       spacing in x and z.
    /// @throws std::invalid_argument if geometry does not have the grid size
    ///         and grid spacing.
    void initialize(const Geometry2D &geometry);
    /// @result True indicates that the class is initialized.
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @name Grid Geometry Information
    /// @{
    /// @result The number of grid points.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfGridPoints() const; 
    /// @result The number of grid points.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] int getNumberOfCells() const;
    /// @result The model geometry.
    /// @throws std::runtime_error if \c isInitialized() is false.
    [[nodiscard]] Geometry2D getGeometry() const;
    /// @}

    /// @name Velocity and Slowness Models
    /// @{
    /// @brief Sets a nodal velocity model.
    /// @param[in] nGrid  The number of grid points where nGrid = nGridZ*nGridX.
    ///                   The number of grid points in each direction can be
    ///                   obtained from \c getGeometry().
    /// @param[in] velocity  The velocity model in m/s defined at the
    ///                      the grid points.  This is an [nGridZ x nGridX] or
    ///                      [nGridX x nGridZ] matrix depending on the ordering.
    /// @param[in] ordering  The ordering of the velocity model.
    /// @throws std::invalid_argument if nGrid != \c getNumberOfGridPoints()
    ///         or velocity is NULL or has a non-positive value.
    /// @throws std::runtime_error if the class is not initialized.
    void setNodalVelocities(int nGrid, const T velocity[],
                            EikonalXX::Ordering2D ordering);
 
    /// @brief Sets the cell-based velocity model.
    /// @param[in] nCell  The number of cells where nCell = nCellZ*nCellX.
    ///                   The number of cells in each direction can be obtained
    ///                   obtained from \c getGeometry().
    /// @param[in] velocity  The velocity model in m/s defined at the cell
    ///                      centers.  This is an [nCellZ x nCellX] or 
    ///                      [nCellX x nCellZ] matrix depending on the ordering.
    /// @param[in] ordering  The ordering of the velocity model.
    /// @throws std::invalid_argument if nCell != \c getNumberOfCells()
    ///         or velocity is NULL or has a non-positive value.
    /// @throws std::runtime_error if the class is not initialized.
    template<typename U>
    void setCellularVelocities(int nCell, const U velocity[],
                               EikonalXX::Ordering2D ordering); 

    /// @result True indicates the velocity model was set.
    [[nodiscard]] bool haveVelocities() const noexcept;

    /// @result The cell-based slowness model in seconds/meter in natural
    ///         ordering.
    /// @throws std::runtime_error if the velocity model was not set.
    /// @sa \c haveVelocities()
    [[nodiscard]] std::vector<T> getSlowness() const;
    /// @result A pointer to the cell-based slowness model in seconds/meter
    ///         in natural ordering.
    /// @throws std::runtime_error if the velocity model was not set.
    /// @sa \c haveVelocities()
    [[nodiscard]] const T* getSlownessPointer() const;
    /// @brief Writes the velocity model to VTK.
    /// @param[in] fileName  The name of the VTK file.
    /// @param[in] title     The dataset's title.
    /// @throws std::runtime_error if \c haveVelocities() is false.
    /// @throws std::invalid_argument if there is an error while opening the
    ///         output file.
    void writeVTK(const std::string &fileName,
                  const std::string &title = "velocity_m/s") const;
  
    /// @result The cell-based velocity model in meters/second in natural
    ///         ordering.
    /// @throws std::runtime_error if \c haveVelocities() is false.
    [[nodiscard]] std::vector<T> getVelocities() const;
    /// @}  
private:
    class Model2DImpl;
    std::unique_ptr<Model2DImpl> pImpl;
};
}
#endif
