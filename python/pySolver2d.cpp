#include "eikonalxx/solver2d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/station2d.hpp"
#include "include/pySolver2d.hpp"
#include "include/pyGeometry2d.hpp"
#include "include/pyModel2d.hpp"
#include "include/pySolverOptions.hpp"
#include "include/pySource2d.hpp"
#include "include/pyStation2d.hpp"
#include "include/pySolverOptions.hpp"

using namespace PEikonalXX;

/// Constructor
Solver2D::Solver2D() :
    mSolver(std::make_unique<EikonalXX::Solver2D<double>> ())
{
}

/// Destructor
Solver2D::~Solver2D() = default;


void Solver2D::initialize(const SolverOptions &options,
                          const Geometry2D &geometry)
{
    mSolver->initialize(*options.getNativeClassPointer(),
                        *geometry.getNativeClassPointer());
}

SolverOptions Solver2D::getOptions() const
{
    return SolverOptions{mSolver->getOptions()};
}

Geometry2D Solver2D::getGeometry() const
{
    return Geometry2D{mSolver->getGeometry()};
}

bool Solver2D::isInitialized() const noexcept
{
    return mSolver->isInitialized();
}

void Solver2D::setVelocityModel(const Model2D &velocityModel)
{
    mSolver->setVelocityModel(*velocityModel.getNativeClassPointer());
}

Model2D Solver2D::getVelocityModel() const
{
    return Model2D{mSolver->getVelocityModel()};
}

bool Solver2D::haveVelocityModel() const noexcept
{
    return mSolver->haveVelocityModel();
}

double Solver2D::getSlowness(const int iCellX, const int iCellZ) const
{
    return static_cast<double> (mSolver->getSlowness(iCellX, iCellZ));
}

void Solver2D::setSource(const Source2D &source)
{
    mSolver->setSource(*source.getNativeClassPointer());
}

Source2D Solver2D::getSource() const
{
    return Source2D {mSolver->getSource()};
}

bool Solver2D::haveSource() const noexcept
{
    return mSolver->haveSource();
}

void Solver2D::solve() 
{
    mSolver->solve();
}

//std::vector<double> getTravelTimeField() const;
bool Solver2D::haveTravelTimeField() const noexcept
{
    return mSolver->haveTravelTimeField();
}

//std::vector<double> getTravelTimeGradientField() const;
bool Solver2D::haveTravelTimeGradientField() const noexcept
{
    return mSolver->haveTravelTimeGradientField();
}



