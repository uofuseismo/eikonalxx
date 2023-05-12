#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
#include <vector>
namespace EikonalXX::AbstractBaseClass
{
template<class T>
class ISolver2D
{
public:
    virtual ~ISolver2D() = default;
    virtual std::vector<T> getTravelTimeField() const = 0;
    virtual const T* getTravelTimeFieldPointer() const = 0;
    virtual std::vector<T> getTravelTimeGradientFieldInX() const = 0;
    virtual const T* getTravelTimeGradientFieldInXPointer() const = 0;
    virtual std::vector<T> getTravelTimeGradientFieldInZ() const = 0;
    virtual const T* getTravelTimeGradientFieldInZPointer() const = 0;
    virtual bool haveTravelTimeField() const noexcept = 0;
    virtual bool haveTravelTimeGradientField() const noexcept = 0;
};
}
#endif
