#ifndef EIKONALXX_IO_VTKRECTILINEARGRID2D_HPP
#define EIKONALXX_IO_VTKRECTILINEARGRID2D_HPP
#include <string>
#include <memory>
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

    /// @brief Sets the 2D geometry.
    void setGeometry(const Geometry2D &geometry);

    /// @brief Adds a nodal dataset.
    /// @param[in] name     The name of the dataset.
    /// @param[in] dataSet  The node-based dataset.  This is stored in
    ///                     EikonalXX's natural ordering.
    void addNodalDataset(const std::string &name, const double *dataSet);
    /// @copydoc addNodalDataset 
    void addNodalDataset(const std::string &name, const float *dataSet);

    /// @brief Adds a cell-based dataset.
    /// @param[in] name     The name of the dataset.
    /// @param[in] dataSet  The node-based dataset.  This is stored in
    ///                     EikonalXX's natural ordering.

    /// @brief Releases references to all datasets.
    void clearDatasets();
    /// @name Destructors
    /// @{
    ~VTKRectilinearGrid2D();
    /// @}
private:
    class VTKRectlinearGridImpl;
    std::unique_ptr<VTKRectlinearGridImpl> pImpl;
};
}
}
#endif
