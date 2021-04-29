#include <cmath>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include <gtest/gtest.h>

namespace
{
TEST(Analytic, Homogeneous2D)
{
    double velocity = 3000;
    int nx = 55;
    int nz = 23;
    double dx = 100;
    double dz = 100;
    double x0 = 1;
    double z0 = 2;
    double xSrc = (nx - 1)/4.*dx + x0 + 4;
    double zSrc = (nz - 1)/3.*dz + z0 - 40;   
    // Create geometry
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
    // Initialize solver
    EikonalXX::Analytic::Homogeneous2D<double> solver;
    EXPECT_NO_THROW(solver.initialize(geometry));
    EXPECT_NO_THROW(solver.setSource(std::pair(xSrc, zSrc)));
    EXPECT_NO_THROW(solver.setVelocityModel(velocity));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_TRUE(solver.haveSource());
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_NO_THROW(solver.solve());
    auto travelTimes = solver.getTravelTimeField();
    //solver.writeVTK("test.vtk");
    // Verify
    double dtMax = 0;
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            auto x = x0 + ix*dx;
            auto z = z0 + iz*dz;
            auto delX = x - xSrc;
            auto delZ = z - zSrc;
            auto t = std::hypot(delX, delZ)/velocity;
            auto indx = iz*nx + ix;
            dtMax = std::max(dtMax, std::abs(t - travelTimes.at(indx)));
        }
    } 
    EXPECT_NEAR(dtMax, 0, 1.e-15);
}     
}
