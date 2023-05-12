#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
namespace EikonalXX::AbstractBaseClass
{
template<class T>
class ISolver3D
{
public:
    virtual ~ISolver3D() = default;
    virtual std::vector<T> getTravelTimeField() const = 0;
    virtual const T* getTravelTimeFieldPointer() const = 0;
    virtual std::vector<T> getTravelTimeGradientFieldInX() const = 0;
    virtual const T* getTravelTimeGradientFieldInXPointer() const = 0;
    virtual std::vector<T> getTravelTimeGradientFieldInY() const = 0;
    virtual const T* getTravelTimeGradientFieldInYPointer() const = 0;
    virtual std::vector<T> getTravelTimeGradientFieldInZ() const = 0;
    virtual const T* getTravelTimeGradientFieldInZPointer() const = 0;
    virtual bool haveTravelTimeField() const noexcept = 0;
    virtual bool haveTravelTimeGradientField() const noexcept = 0;
};
}
#endif
