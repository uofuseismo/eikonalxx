#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "gradient.hpp"
#include <gtest/gtest.h>

namespace
{


TEST(Gradient, Homogeneous2D)
{
    EikonalXX::Analytic::Homogeneous2D<double> solver;
    const double velocity{2500};
    const int nx = 201;
    const int nz = 51;
    double dx = 10;
    double dz = 10.5;
    double x0 = 0;
    double z0 = 0;
    double xSrc = (nx - 1)/4.*dx + x0 + 0.3*dx;
    double zSrc = (nz - 1)/3.*dz + z0 - 0.4*dx;
    // Create geometry
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
    // Initialize solver
    solver.initialize(geometry);
    solver.setVelocityModel(velocity);  
    solver.setSource(std::pair {xSrc, zSrc});
    // Solve eikonal equation and compute analytic travel time fields
    solver.solve();
    solver.computeTravelTimeGradientField();
    // Pass this to the gradient
}

}
