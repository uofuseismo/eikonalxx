#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
namespace EikonalXX
{
 class Source3D;
}
namespace EikonalXX::AbstractBaseClass
{
template<class T>
/// @class ISolver3D "solver3d.hpp" "eikonalxx/abstractBaseClass/solver3d.hpp"
/// @brief The abstract base class for the three dimensional solvers.
/// @ingroup Solver_BaseClass
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class ISolver3D
{
public:
    /// @brief Default destructor.
    virtual ~ISolver3D() = default;
    /// @result True indicates the velocity model was set.
    virtual bool haveVelocityModel() const noexcept = 0;
    /// @result The travel time field.
    virtual std::vector<T> getTravelTimeField() const = 0;
    /// @result A pointer to the travel time field.
    virtual const T* getTravelTimeFieldPointer() const = 0;
    /// @result The gradient of the travel time field.
    virtual std::vector<T> getTravelTimeGradientField() const = 0;
    /// @result The pointer to the gradient of the travel time field.
    virtual const T* getTravelTimeGradientFieldPointer() const = 0;
    /// @result True indicates the travel time field is available.
    virtual bool haveTravelTimeField() const noexcept = 0;
    /// @result True indicates the gradient of the travel time field
    ///         is available.
    virtual bool haveTravelTimeGradientField() const noexcept = 0;
    /// @result The source information.
    virtual Source3D getSource() const = 0;
    /// @result True indicates the source exists.
    virtual bool haveSource() const noexcept = 0;
    /// @result The slowness at a cell.
    virtual T getSlowness(int iCellX, int iCellY, int iCellZ) const = 0;
};
}
#endif
