#ifndef EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
#define EIKONALX_ABSTRACT_BASE_CLASS_SOLVER_2D_HPP
namespace EikonalXX::AbstractBaseClass
{
template<class T>
class ISolver2D
{
public:
    virtual std::vector<T> getTravelTimeField() const = 0;
    virtual const T* getTravelTimeFieldPointer() const = 0;
    virtual bool haveTravelTimeField() const noexcept = 0;
};
}
#endif
