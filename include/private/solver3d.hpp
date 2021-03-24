#ifndef EIKONALXX_PRIVATE_SOLVER3D_HPP
#define EIKONALXX_PRIVATE_SOLVER3D_HPP
#define HUGE 1.e10
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/graph3d.hpp"
#include "solverUtilities3d.hpp"

namespace
{

using namespace EikonalXX;

/// Solver class
template<class T, EikonalXX::SweepNumber3D E>
class Solver3DSweep
{
public:
    /// Initialize the variables in this class 
    void initialize(const EikonalXX::Geometry3D &geometry,
                    const EikonalXX::SolverOptions &options)
    {
        mX = geometry.getNumberOfGridPointsInX();
        mY = geometry.getNumberOfGridPointsInY();
        mZ = geometry.getNumberOfGridPointsInZ();
        mGraph.initialize(mX, mY, mZ);
    }
    /// Clears the class / releases memory
    void clear() noexcept
    {
        mGraph.clear();
    }
//private:
    Graph3D<E> mGraph;
    /// Number of grid points in x.
    int mX = 0;
    int mY = 0;
    int mZ = 0;
};

///--------------------------------------------------------------------------///
///                               Template Instantiation                     ///
///--------------------------------------------------------------------------///

template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP1>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP2>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP3>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP4>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP5>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP6>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP7>;
template class Solver3DSweep<double, EikonalXX::SweepNumber3D::SWEEP8>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP1>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP2>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP3>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP4>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP5>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP6>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP7>;
template class Solver3DSweep<float, EikonalXX::SweepNumber3D::SWEEP8>;


}
#endif
