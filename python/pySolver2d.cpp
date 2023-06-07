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

void Solver2D::writeVTK(const std::string &fileName,
                        const std::string &title,
                        const bool writeGradientField) const
{
    mSolver->writeVTK(fileName, title, writeGradientField);
}

pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> 
Solver2D::getTravelTimeField() const
{
    const auto *travelTimes = mSolver->getTravelTimeFieldPointer();
    auto nGrid = mSolver->getGeometry().getNumberOfGridPoints();
    auto result = pybind11::array_t<double, pybind11::array::c_style> (nGrid);
    pybind11::buffer_info resultBuffer = result.request();
    auto resultPointer = static_cast<double *> (resultBuffer.ptr);
    std::copy(travelTimes, travelTimes + nGrid, resultPointer);
    return result;
}

pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast> 
Solver2D::getTravelTimeGradientField() const
{
    const auto *gradientTravelTimes
        = mSolver->getTravelTimeGradientFieldPointer();
    auto nGrid2 = 2*mSolver->getGeometry().getNumberOfGridPoints();
    auto result = pybind11::array_t<double, pybind11::array::c_style> (nGrid2);
    pybind11::buffer_info resultBuffer = result.request();
    auto resultPointer = static_cast<double *> (resultBuffer.ptr);
    std::copy(gradientTravelTimes, gradientTravelTimes + nGrid2, resultPointer);
    return result;
}

void Solver2D::computeTravelTimeGradientField()
{
    mSolver->computeTravelTimeGradientField();
}

bool Solver2D::haveTravelTimeField() const noexcept
{
    return mSolver->haveTravelTimeField();
}

bool Solver2D::haveTravelTimeGradientField() const noexcept
{
    return mSolver->haveTravelTimeGradientField();
}

void Solver2D::setStations(const std::vector<Station2D> &stationsIn)
{
    if (stationsIn.empty()){return;}
    std::vector<EikonalXX::Station2D> stations;
    stations.reserve(stationsIn.size());
    for (const auto &station : stationsIn)
    {
        stations.push_back(*station.getNativeClassPointer());
    }
    mSolver->setStations(stations);
}

std::vector<Station2D> Solver2D::getStations() const
{
    auto stations = mSolver->getStationsReference();
    std::vector<Station2D> result;
    for (const auto &station : stations)
    {
        result.push_back(station);
    }
    return result;
}

/// Initialize class
void PEikonalXX::initializeSolver2D(pybind11::module &m) 
{
    pybind11::class_<PEikonalXX::Solver2D> s(m, "Solver2D");
    s.def(pybind11::init<> ());
    s.doc() = R"""(
This defines the 2D solver.  The usage is as follows:

Step 1: Initialize the solver with the initialize method.
Step 2: Set the velocity model with the velocity_model property.
Step 3: Set the source with the source property.
Step 4: Optionally, set the station positions with the station property.
Step 5: Solve for the eikonal equation with the solve method.
Step 6: Optionally compute the gradient of the travel time field with
        the compute_travel_time_gradient_field method.

The solver was designed so that initializiation need only be performed once. 
After that, you can change the velocity model, source, and stations as much
as you like.

Properties
----------
source : Source2D
    Defines the seismic source in the model.  This must be called after initialization.
velocity_model : VelocityModel2D
    Defines the velocity model.  This must be called after initialization.

Read-Only Properties
--------------------
geometry : Geometry2D
    The 2D geometry.
have_source : bool
    True indicates the source was set.
have_travel_time_field : bool
    True indicates the travel time field was computed.
have_velocity_model : bool
    True indicates the velocity model was set.
travel_time_field : np.array
    The travel time field in seconds.  You can reshape this with:
    result.reshape(geometry.nz, geometry.nx)
travel_time_gradient_field : np.array
    The gradient of the travel time field in seconds/meter in x and z.
    You can reshape this with: result.reshape(geometry.nz, geometry.nx, 2)

)""";
    s.def("initialize",
          &Solver2D::initialize,
          "Initializes the solver.");
    s.def("solve",
          &Solver2D::solve,
          "Solves the eikonal equation.  Note, the class must be initialized and the source and velocity model must be set.");
    s.def("compute_travel_time_gradient_field",
          &Solver2D::computeTravelTimeGradientField,
          "Computes the gradient of the travel time field.");
    s.def("write_vtk",
          &Solver2D::writeVTK,
          "Writes the travel time field and, optionally, the gradient field to disk.  Note, have_travel_time_field must be true.",
          pybind11::arg("file_name"),
          pybind11::arg("title") = "travel_time_field_s",
          pybind11::arg("write_gradient_field") = false);
    s.def_property("source",
                   &Solver2D::getSource,
                   &Solver2D::setSource);
    s.def_property("velocity_model",
                   &Solver2D::getVelocityModel,
                   &Solver2D::setVelocityModel);
    s.def_property("stations",
                   &Solver2D::getStations,
                   &Solver2D::setStations);
    s.def_property_readonly("have_source",
                            &Solver2D::haveSource);
    s.def_property_readonly("have_velocity_model",
                            &Solver2D::haveVelocityModel);
    s.def_property_readonly("have_travel_time_field",
                            &Solver2D::haveTravelTimeField);
    s.def_property_readonly("have_travel_time_gradient_field",
                            &Solver2D::haveTravelTimeGradientField);
    s.def_property_readonly("travel_time_field",
                            &Solver2D::getTravelTimeField);
    s.def_property_readonly("travel_time_gradient_field",
                            &Solver2D::getTravelTimeGradientField);
}

