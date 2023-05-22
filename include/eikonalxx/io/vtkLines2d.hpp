#ifndef EIKONALXX_IO_VTK_LINES_2D_HPP
#define EIKONALXX_IO_VTK_LINES_2D_HPP
#include <string>
#include <vector>
#include <memory>
#include "eikonalxx/enums.hpp"
namespace EikonalXX
{
class Geometry2D;
}
namespace EikonalXX::IO
{
/// @class VTKLines2D "vtkLines2d.hpp" "eikonalxx/io/vtkLines2d.hpp"
/// @brief A utility class for writing polygons to Legacy VTK file format.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class VTKLines2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    VTKLines2D();
    /// @}

    /// @name Step 1: Open the VTK file 
    /// @{

    /// @brief Opens a VTK file for writing.
    /// @param[in] fileName  The name of the VTK file to write to.
    /// @param[in] geometry  The model geometry.
    /// @param[in] title     The title for the dataset.
    /// @throws std::invalid_argument if the file path cannot be created.
    void open(const std::string &fileName,
              const EikonalXX::Geometry2D &geometry,
              const std::string &title = "ray_paths");
    /// @result True indicates that the the file is open for writing.
    [[nodiscard]] bool isOpen() const noexcept;
    /// @}

    /// @name Step 2: Write datasets
    /// @{

    /// @brief Writes a nodal dataset.
    /// @param[in] polygons  The polygons to write.
    /// @throws std::runtime_error if \c isOpen() is false.
    void write(const std::vector<std::vector< std::pair<double, double>> > &polygons) const;
    /// @}

    /// @name Step 3: Close the file
    /// @{

    /// @brief Closes the file.
    void close() noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~VTKLines2D();
    /// @}

    // Remove some functionality
    VTKLines2D(const VTKLines2D &vtk) = delete;
    VTKLines2D(VTKLines2D &&vtk) noexcept = delete;
    VTKLines2D& operator=(const VTKLines2D &vtk) = delete;
    VTKLines2D& operator=(VTKLines2D &&vtk) noexcept = delete;
private:
    class VTKLines2DImpl;
    std::unique_ptr<VTKLines2DImpl> pImpl;
};
}
#endif
