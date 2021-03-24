#include "eikonalxx/solver3d.hpp"
#include "eikonalxx/graph3d.hpp"
#include "private/solver3d.hpp"

using namespace EikonalXX;

template<class T>
class Solver3D<T>::Solver3DImpl
{
public:
};

/// C'tor 
template<class T>
Solver3D<T>::Solver3D() :
    pImpl(std::make_unique<Solver3DImpl> ())
{
    Solver3DSweep<T, SweepNumber3D::SWEEP1> mSweep1;
    Solver3DSweep<T, SweepNumber3D::SWEEP1> mSweep2;
}

/// Destructor
template<class T>
Solver3D<T>::~Solver3D() = default;

///--------------------------------------------------------------------------///
///                           Template Class Instantiation                   ///
///--------------------------------------------------------------------------///
template class EikonalXX::Solver3D<double>;
template class EikonalXX::Solver3D<float>;
