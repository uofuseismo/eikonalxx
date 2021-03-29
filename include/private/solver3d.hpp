#ifndef EIKONALXX_PRIVATE_SOLVER3D_HPP
#define EIKONALXX_PRIVATE_SOLVER3D_HPP
#include <limits>
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/graph3d.hpp"
#include "solverUtilities3d.hpp"


#define HUGE 1.e10

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
        // Uniform grid?
        auto dx = geometry.getGridSpacingInX();
        auto dy = geometry.getGridSpacingInY();
        auto dz = geometry.getGridSpacingInZ();
        mUniformGrid = false;
        if (std::abs(dx - dy) < std::numeric_limits<T>::epsilon()*100 &&
            std::abs(dx - dz) < std::numeric_limits<T>::epsilon()*100)
        {
            mUniformGrid = true;
        }
        mDx = static_cast<T> (dx);
        mDy = static_cast<T> (dy);
        mDz = static_cast<T> (dz);
        // Initialize the graph
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
    /// Grid spacing in x, y, and z
    float mDx = 0;
    float mDy = 0;
    float mDz = 0;
    /// Number of grid points in x, y, and z.
    int mX = 0;
    int mY = 0;
    int mZ = 0;
    bool mUniformGrid = true;
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
