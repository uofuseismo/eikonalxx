#ifndef EIKONALXX_IO_VTKRECTILINEARGRID3D_HPP
#define EIKONALXX_IO_VTKRECTILINEARGRID3D_HPP
#include <string>
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
class Geometry3D;
namespace IO
{
class VTKRectilinearGrid3D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    VTKRectilinearGrid3D();
    /// @}

    /// @name Destructors
    /// @{
    ~VTKRectilinearGrid3D();
    /// @}
    
    /// @name Step 1: Open the VTK file 
    /// @{
    /// @brief Opens a VTK file for writing.
    /// @param[in] fileName     The name of the VTK file to write to.
    /// @param[in] geometry     The 3D geometry.
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
              const Geometry3D &geometry,
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
                           Ordering3D ordering = Ordering3D::NATURAL) const;
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
                              Ordering3D ordering = Ordering3D::NATURAL) const;
    /// @}

    /// @name Step 3: Close the file
    /// @{
    /// @brief Closes the file.
    void close() noexcept;
    /// @}

    // Remove some functionality
    VTKRectilinearGrid3D(const VTKRectilinearGrid3D &vtk) = delete;
    VTKRectilinearGrid3D(VTKRectilinearGrid3D &&vtk) noexcept = delete;
    VTKRectilinearGrid3D& operator=(const VTKRectilinearGrid3D &vtk) = delete;
    VTKRectilinearGrid3D& operator=(VTKRectilinearGrid3D &&vtk) noexcept = delete;
private:
    class VTKRectilinearGrid3DImpl;
    std::unique_ptr<VTKRectilinearGrid3DImpl> pImpl;
};
}
}
#endif
