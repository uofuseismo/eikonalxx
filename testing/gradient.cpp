#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/analytic/linearGradient2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/model2d.hpp"
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
    // Create velocity model
    EikonalXX::Model2D<double> model;
    std::vector<double> vModel((nx - 1)*(nz - 1), velocity); 
    model.initialize(geometry);
    model.setCellularVelocities(vModel.size(), vModel.data(),
                                EikonalXX::Ordering2D::Natural);
    // Initialize solver
    solver.initialize(geometry);
    solver.setVelocityModel(velocity);  
    solver.setSource(std::pair {xSrc, zSrc});
    // Solve eikonal equation and compute analytic travel time fields
    solver.solve();
    solver.computeTravelTimeGradientField();
    //solver.writeVTK("debug.vtk", "homogeneous_field", true);
    auto gradReference = solver.getTravelTimeGradientField();
    // Pass this to the gradient
    std::vector<double> gradient;
    ::finiteDifference(geometry, solver.getSource(), model,
                       solver.getTravelTimeFieldPointer(),
                       &gradient,
                       ::DerivativeType::CentralDifference);
//                      const DerivativeType derivativeType = DerivativeType::CentralDifference)
    EXPECT_EQ(gradient.size(), gradReference.size());
    for (int i = 0; i < nx*nz; ++i)
    {
        double resx = gradReference[2*i]     - gradient[2*i]; 
        double resz = gradReference[2*i + 1] - gradient[2*i + 1];
        EXPECT_NEAR(std::abs(resx), 0, 5.e-5);
        EXPECT_NEAR(std::abs(resz), 0, 5.e-5);
    }
    /*
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            auto index = 2*::gridToIndex(nx, ix, iz);
            double res = gradReference[index] - gradient[index];
            if (std::abs(res) > 1.e-4)
            {   
                std::cout << ix << " " << iz << " " << gradient[index] << " " << gradReference[index] << std::endl;
            }
        }
    }
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            auto index = 2*::gridToIndex(nx, ix, iz) + 1;
            double res = gradReference[index] - gradient[index];
            if (std::abs(res) > 1.e-4)
            {
                std::cout << ix << " " << iz << " " << gradient[index] << " " << gradReference[index] << std::endl;
            }
        }
    }
    */
}

TEST(Gradient, LinearGradient2D)
{
    std::pair<double, double> velocity({1000, 4000});
    int nx = 401; 
    int nz = 201; 
    double dx = 100;
    double dz = 101;
    double x0 = 0;
    double z0 = 0;
    double xSrc = (nx - 1)/2.*dx + x0;
    double zSrc = z0;
    // Create geometry
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
    // Compute the model
    EikonalXX::Model2D<double> model;
    std::vector<double> vModel(nx*nz, 0);
    for (int iz = 0; iz < nz; ++iz)
    {
        double zi = z0 + iz*dz;
        auto vi = velocity.first
                + (zi - z0)*(velocity.second - velocity.first)
                 /( (nz - 1)*dz );
        for (int ix = 0; ix < nx; ++ix)
        {
            vModel.at(::gridToIndex(nx, ix, iz)) = vi; 
        }
    }
    model.initialize(geometry);
    model.setNodalVelocities(vModel.size(), vModel.data(),
                             EikonalXX::Ordering2D::Natural);
    // Initialize solver
    EikonalXX::Analytic::LinearGradient2D<double> solver;
    EXPECT_NO_THROW(solver.initialize(geometry));
    EXPECT_NO_THROW(solver.setSource(std::pair(xSrc, zSrc)));
    EXPECT_NO_THROW(solver.setVelocityModel(velocity));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_TRUE(solver.haveSource());
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_NO_THROW(solver.solve());
    solver.computeTravelTimeGradientField();
    //solver.writeVTK("linearGradient.vtk", "linear_gradient_field", true);
    auto gradReference = solver.getTravelTimeGradientField();
    // Pass this to the gradient
    std::vector<double> gradient;
    ::finiteDifference(geometry, solver.getSource(), model,
                       solver.getTravelTimeFieldPointer(),
                       &gradient,
                       ::DerivativeType::CentralDifference);
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            int i = 2*::gridToIndex(nx, ix, iz);
            double resx = gradReference[i]     - gradient[i];
            double resz = gradReference[i + 1] - gradient[i + 1];
            EXPECT_NEAR(std::abs(resx), 0, 1.e-4);
            // Near the source the finite difference is results in a 
            // positive number.  This is because the travel time at the grid
            // point diagonal the source is bigger than the travel time at
            // at the adjacent grid point (e.g., du/dz ~ u_{z+1} - u_z).
            // However, the analytic method actually can account for the fact
            // that the ray is actually curving upward in the cell!  Hence,
            // why we don't compare these two points - basically the analytic
            // solution is too good.
            if (!( (iz == 0 && ix == 199) ||
                   (iz == 0 && ix == 200) ||
                   (iz == 0 && ix == 201)) )
            {
                EXPECT_NEAR(std::abs(resz), 0, 1.e-4);
            }
        }
    }
}

}
