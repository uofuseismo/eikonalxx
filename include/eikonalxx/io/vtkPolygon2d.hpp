#ifndef EIKONALXX_IO_VTK_POLYGON_2D_HPP
#define EIKONALXX_IO_VTK_POLYGON_2D_HPP
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
/// @class VTKPolygon2D "vtkPolygon2d.hpp" "eikonalxx/io/vtkPolygon2d.hpp"
/// @brief A utility class for writing polygons to Legacy VTK file format.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class VTKPolygon2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    VTKPolygon2D();
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
    ~VTKPolygon2D();
    /// @}

    // Remove some functionality
    VTKPolygon2D(const VTKPolygon2D &vtk) = delete;
    VTKPolygon2D(VTKPolygon2D &&vtk) noexcept = delete;
    VTKPolygon2D& operator=(const VTKPolygon2D &vtk) = delete;
    VTKPolygon2D& operator=(VTKPolygon2D &&vtk) noexcept = delete;
private:
    class VTKPolygon2DImpl;
    std::unique_ptr<VTKPolygon2DImpl> pImpl;
};
}
#endif
