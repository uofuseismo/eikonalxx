#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_3D_HPP
namespace EikonalXX::AbstractBaseClass
{
template<class T>
class ISolver3D
{
public:
    virtual std::vector<T> getTravelTimeField() const = 0;
    virtual const T* getTravelTimeFieldPointer() const = 0;
    virtual bool haveTravelTimeField() const noexcept = 0;
};
}
#endif
