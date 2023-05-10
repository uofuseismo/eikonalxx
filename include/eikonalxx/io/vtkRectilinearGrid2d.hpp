#ifndef EIKONALXX_IO_VTKRECTILINEARGRID2D_HPP
#define EIKONALXX_IO_VTKRECTILINEARGRID2D_HPP
#include <string>
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
class Geometry2D;
namespace IO
{
class VTKRectilinearGrid2D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    VTKRectilinearGrid2D();
    /// @}

    /// @name Destructors
    /// @{
    ~VTKRectilinearGrid2D();
    /// @}
    
    /// @name Step 1: Open the VTK file 
    /// @{
    /// @brief Opens a VTK file for writing.
    /// @param[in] fileName     The name of the VTK file to write to.
    /// @param[in] geometry     The 2D geometry.
    /// @param[in] title        Specifies the title.  If this exceeds 256
    ///                         characters then it will be truncated.  
    ///                         Additionally, blank spaces will be filled by
    ///                         underscores. 
    /// @param[in] writeBinary  If true then the output will be binary instead
    ///                         of plain text.  
    /// @throws std::invalid_argument if the number of grid points or grid
    ///         spacing is not fully specified on the geometry or if the file
    ///         path cannot be created. 
    void open(const std::string &fileName,
              const Geometry2D &geometry,
              const std::string &title = "traveltimes",
              const bool writeBinary = true);
    /// @result True indicates that the the file is open for writing.
    [[nodiscard]] bool isOpen() const noexcept;
    /// @}

    /// @name Step 2: Write datasets
    /// @{
    /// @brief Writes a nodal dataset.
    /// @param[in] name      The name of the nodal dataset.  Note, all blank
    ///                      spaces will be replaced with underscores.
    /// @param[in] data      The nodal dataset to write.  This has dimension
    ///                      geometry.getNumberOfGridPoints().
    /// @param[in] ordering  The dataset's ordering.
    /// @throws std::runtime_error if \c isOpen() is false.
    /// @throws std::invalid_argument if the dataset's name is empty or the
    ///         data pointer is NULL.
    template<typename T>
    void writeNodalDataset(const std::string &name,
                           const T *data,
                           Ordering2D ordering = Ordering2D::Natural) const;
    /// @brief Writes a cell-based dataset.
    /// @brief Writes a nodal dataset.
    /// @param[in] name      The name of the cellular dataset.  Note, all blank
    ///                      spaces will be replaced with underscores.
    /// @param[in] data      The cellular dataset to write.  This has dimension
    ///                      geometry.getNumberOfCells().
    /// @param[in] ordering  The dataset's ordering.
    /// @throws std::runtime_error if \c isOpen() is false.
    /// @throws std::invalid_argument if the dataset's name is empty or the
    ///         data pointer is NULL.
    template<typename T>
    void writeCellularDataset(const std::string &fname,
                              const T *data,
                              Ordering2D ordering = Ordering2D::Natural) const;
    /// @}

    /// @name Step 3: Close the file
    /// @{
    /// @brief Closes the file.
    void close() noexcept;
    /// @}

    // Remove some functionality
    VTKRectilinearGrid2D(const VTKRectilinearGrid2D &vtk) = delete;
    VTKRectilinearGrid2D(VTKRectilinearGrid2D &&vtk) noexcept = delete;
    VTKRectilinearGrid2D& operator=(const VTKRectilinearGrid2D &vtk) = delete;
    VTKRectilinearGrid2D& operator=(VTKRectilinearGrid2D &&vtk) noexcept = delete;
private:
    class VTKRectilinearGrid2DImpl;
    std::unique_ptr<VTKRectilinearGrid2DImpl> pImpl;
};
}
}
#endif
