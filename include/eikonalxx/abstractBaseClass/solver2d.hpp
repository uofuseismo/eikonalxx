#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
#include <vector>
namespace EikonalXX
{
 class Source2D;
}
namespace EikonalXX::AbstractBaseClass
{
template<class T>
/// @class ISolver2D "solver2d.hpp" "eikonalxx/abstractBaseClass/solver2d.hpp"
/// @brief The abstract base class for the two dimensional solvers.
/// @ingroup Solver_BaseClass
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class ISolver2D
{
public:
    /// @brief Default destructor.
    virtual ~ISolver2D() = default;
    /// @result True indicates the velocity model was set.
    virtual bool haveVelocityModel() const noexcept = 0;
    /// @brief Gets a copy of the travel time field.
    virtual std::vector<T> getTravelTimeField() const = 0;
    /// @brief Gets a pointer to the travel time field.
    virtual const T* getTravelTimeFieldPointer() const = 0;
    /// @brief Gets a copy of the gradient of the travel time field.
    virtual std::vector<T> getTravelTimeGradientField() const = 0;
    /// @brief Gets a pointer to the gradient of the travel time field.
    virtual const T *getTravelTimeGradientFieldPointer() const = 0;
    /// @result True indicates the travel time field was computed.
    virtual bool haveTravelTimeField() const noexcept = 0;
    /// @result True indicates the gradient of the travel time field was computed. 
    virtual bool haveTravelTimeGradientField() const noexcept = 0;
    /// @result The source information.
    virtual Source2D getSource() const = 0;
    /// @result True indicates the source exists.
    virtual bool haveSource() const noexcept = 0;
    /// @result The slowness at a cell.
    virtual T getSlowness(int iCellX, int iCellZ) const = 0;
};
}
#endif
