#include "eikonalxx/solver3d.hpp"
#include "eikonalxx/source3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/graph3d.hpp"
#include "eikonalxx/model3d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include "private/solver3d.hpp"

using namespace EikonalXX;

template<class T>
class Solver3D<T>::Solver3DImpl
{
public:
    Solver3DSweep<T, SweepNumber3D::SWEEP1> mSolverSweep1;
    Solver3DSweep<T, SweepNumber3D::SWEEP2> mSolverSweep2;
    Solver3DSweep<T, SweepNumber3D::SWEEP3> mSolverSweep3;
    Solver3DSweep<T, SweepNumber3D::SWEEP4> mSolverSweep4;
    Solver3DSweep<T, SweepNumber3D::SWEEP5> mSolverSweep5;
    Solver3DSweep<T, SweepNumber3D::SWEEP6> mSolverSweep6;
    Solver3DSweep<T, SweepNumber3D::SWEEP7> mSolverSweep7;
    Solver3DSweep<T, SweepNumber3D::SWEEP8> mSolverSweep8;
    Model3D<T> mVelocityModel;
    Geometry3D mGeometry;
    Source3D mSource;
    SolverOptions mOptions;
    int mSourceCell = 0;
    bool mHaveTravelTimeField = false;
    bool mHaveVelocityModel = false;
    bool mHaveSource = false;
    bool mInitialized = false;
};

/// C'tor 
template<class T>
Solver3D<T>::Solver3D() :
    pImpl(std::make_unique<Solver3DImpl> ())
{
}

/// Destructor
template<class T>
Solver3D<T>::~Solver3D() = default;

/// Release memory
template<class T>
void Solver3D<T>::clear() noexcept
{
    pImpl->mSolverSweep1.clear();
    pImpl->mSolverSweep2.clear();
    pImpl->mSolverSweep3.clear();
    pImpl->mSolverSweep4.clear();
    pImpl->mVelocityModel.clear();
    pImpl->mGeometry.clear();
    pImpl->mOptions.clear();
    pImpl->mSourceCell = 0;
    pImpl->mHaveTravelTimeField = false;
    pImpl->mHaveVelocityModel = false;
    pImpl->mHaveSource = false;
    pImpl->mInitialized = false;
}

/// Initialize the solver
template<class T>
void Solver3D<T>::initialize(const SolverOptions &options,
                             const Geometry3D &geometry)
{
    clear();
    if (!geometry.haveNumberOfGridPointsInX())
    {
        throw std::invalid_argument("Grid points in x not set on geometry");
    }
    if (!geometry.haveNumberOfGridPointsInY())
    {
        throw std::invalid_argument("Grid points in y not set on geometry");
    }
    if (!geometry.haveNumberOfGridPointsInZ())
    {
        throw std::invalid_argument("Grid points in z not set on geometry");
    }
    if (!geometry.haveGridSpacingInX())
    {
        throw std::invalid_argument("Grid spacing in x not set on geometry");
    }
    if (!geometry.haveGridSpacingInY())
    {
        throw std::invalid_argument("Grid spacing in y not set on geometry");
    }
    if (!geometry.haveGridSpacingInZ())
    {
        throw std::invalid_argument("Grid spacing in z not set on geometry");
    }
    pImpl->mGeometry = geometry;
    pImpl->mOptions = options;
    // Initialize each solver sweep
    pImpl->mSolverSweep1.initialize(options, geometry);
    pImpl->mSolverSweep2.initialize(options, geometry);
    pImpl->mSolverSweep3.initialize(options, geometry);
    pImpl->mSolverSweep4.initialize(options, geometry);
    pImpl->mInitialized = true;
}

/// Initialized?
template<class T>
bool Solver3D<T>::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Sets the velocity model
template<class T>
void Solver3D<T>::setVelocityModel(const Model3D<T> &velocityModel)
{
    pImpl->mHaveVelocityModel = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized())
    {
        throw std::runtime_error("Class not initialized");
    }
    if (!velocityModel.isInitialized())
    {
        throw std::invalid_argument("Velocity model not initialized");
    }
    if (velocityModel.getGeometry() != pImpl->mGeometry)
    {
        throw std::invalid_argument(
            "Velocity model's geometry does not match the solver's geometry");
    }
    if (pImpl->mOptions.getVerbosity() >= Verbosity::INFO)
    {
        std::cout << "Setting velocity model..." << std::endl;
    }
    // Copy the velocity model
    pImpl->mVelocityModel = velocityModel;
    // Set the velocity models for each sweep
    pImpl->mSolverSweep1.setVelocityModel(velocityModel);
    pImpl->mSolverSweep2.setVelocityModel(velocityModel);
    pImpl->mSolverSweep3.setVelocityModel(velocityModel);
    pImpl->mSolverSweep4.setVelocityModel(velocityModel);
/*
    pImpl->mSolverSweep5.setVelocityModel(velocityModel);
    pImpl->mSolverSweep6.setVelocityModel(velocityModel);
    pImpl->mSolverSweep7.setVelocityModel(velocityModel);
    pImpl->mSolverSweep8.setVelocityModel(velocityModel);
*/
    pImpl->mHaveVelocityModel = true;
}

/// Gets the velocity model
template<class T>
Model3D<T> Solver3D<T>::getVelocityModel() const
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    return pImpl->mVelocityModel;
}

/// Have the velocity model?
template<class T>
bool Solver3D<T>::haveVelocityModel() const noexcept
{
    return pImpl->mHaveVelocityModel;
}

/// Set the source from a tuple
template<class T>
void Solver3D<T>::setSource(
    const std::tuple<double, double, double> &sourceLocation)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    auto dx = pImpl->mGeometry.getGridSpacingInX();
    auto dy = pImpl->mGeometry.getGridSpacingInY();
    auto dz = pImpl->mGeometry.getGridSpacingInZ();
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    auto ny = pImpl->mGeometry.getNumberOfGridPointsInY();
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
    auto xmin = pImpl->mGeometry.getOriginInX();
    auto ymin = pImpl->mGeometry.getOriginInY();
    auto zmin = pImpl->mGeometry.getOriginInZ();
    auto xmax = xmin + (nx - 1)*dx;
    auto ymax = ymin + (ny - 1)*dy;
    auto zmax = zmin + (nz - 1)*dz;
    if (std::get<0> (sourceLocation) < xmin ||
        std::get<0> (sourceLocation) > xmax)
    {
        throw std::invalid_argument("x source position = "
                                  + std::to_string(std::get<0> (sourceLocation))
                                  + " must be in range [" + std::to_string(xmin)
                                  + "," + std::to_string(xmax) + "]");
    }
    if (std::get<1> (sourceLocation) < ymin ||
        std::get<1> (sourceLocation) > ymax)
    {
        throw std::invalid_argument("y source position = "
                                  + std::to_string(std::get<1> (sourceLocation))
                                  + " must be in range [" + std::to_string(ymin)
                                  + "," + std::to_string(ymax) + "]");
    }
    if (std::get<2> (sourceLocation) < zmin ||
        std::get<2> (sourceLocation) > zmax)
    {
        throw std::invalid_argument("z source position = "
                                  + std::to_string(std::get<2> (sourceLocation))
                                  + " must be in range [" + std::to_string(zmin)
                                  + "," + std::to_string(zmax) + "]");
    }
    // Tell user what is about to happen
    if (pImpl->mOptions.getVerbosity() >= Verbosity::INFO)
    {
        std::cout << "Setting source location (x,y,z)=("
                  << std::get<0> (sourceLocation) << ","
                  << std::get<1> (sourceLocation) << ","
                  << std::get<2> (sourceLocation) << ")" << std::endl;
    }
    // Get some solver information 
    Source3D source;
    source.setGeometry(pImpl->mGeometry);
    source.setLocationInX(std::get<0> (sourceLocation));
    source.setLocationInY(std::get<1> (sourceLocation));
    source.setLocationInZ(std::get<2> (sourceLocation));
    setSource(source);
}

/// Set the source
template<class T>
void Solver3D<T>::setSource(const Source3D &source)
{
    pImpl->mHaveSource = false;
    pImpl->mHaveTravelTimeField = false;
    if (!isInitialized()){throw std::runtime_error("Class not initialized");}
    if (!source.haveLocationInX())
    {
        throw std::invalid_argument("Source location in x not set");
    }
    if (!source.haveLocationInY())
    {
        throw std::invalid_argument("Source location in y not set");
    }
    if (!source.haveLocationInZ())
    {   
        throw std::invalid_argument("Source location in z not set");
    }   
    pImpl->mSourceCell = source.getCell();
    pImpl->mSource = source;
    pImpl->mHaveSource = true;
    if (pImpl->mOptions.getVerbosity() == Verbosity::DEBUG)
    {
        std::cout << pImpl->mSource << std::endl;
    }
}


/// Get source
template<class T>
Source3D Solver3D<T>::getSource() const
{
    if (!haveSource())
    {
        throw std::runtime_error("Source location not set");
    }
    return pImpl->mSource;
}

/// Have source location?
template<class T>
bool Solver3D<T>::haveSource() const noexcept
{
    return pImpl->mHaveSource;
}

///--------------------------------------------------------------------------///
///                           Template Class Instantiation                   ///
///--------------------------------------------------------------------------///
template class EikonalXX::Solver3D<double>;
template class EikonalXX::Solver3D<float>;
