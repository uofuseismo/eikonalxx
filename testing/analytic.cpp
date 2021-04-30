#include <cmath>
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/analytic/homogeneous3d.hpp"
#include "eikonalxx/analytic/linearGradient2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/source3d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include <gtest/gtest.h>

namespace
{
TEST(Analytic, Homogeneous2D)
{
    double velocity = 3000;
    int nx = 55;
    int nz = 23;
    double dx = 100;
    double dz = 101;
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

TEST(Analytic, Homogeneous3D)
{
    double velocity = 4000;
    int nx = 43; 
    int ny = 43;
    int nz = 23; 
    double dx = 100;
    double dy = 101;
    double dz = 102;
    double x0 = 1;
    double y0 = 1.5;
    double z0 = 2;
    double xSrc = (nx - 1)/4.*dx + x0 + 4;
    double ySrc = 3*(ny - 1)/4.*dy + y0 - 80;
    double zSrc = (nz - 1)/3.*dz + z0 - 40;   
    // Create geometry
    EikonalXX::Geometry3D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInY(ny);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInY(dy);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInY(y0);
    geometry.setOriginInZ(z0);
    // Initialize solver
    EikonalXX::Analytic::Homogeneous3D<double> solver;
    EXPECT_NO_THROW(solver.initialize(geometry));
    EXPECT_NO_THROW(solver.setSource(std::tuple(xSrc, ySrc, zSrc)));
    EXPECT_NO_THROW(solver.setVelocityModel(velocity));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_TRUE(solver.haveSource());
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_NO_THROW(solver.solve());
    auto travelTimes = solver.getTravelTimeField();
    //solver.writeVTK("test3d.vtk");
    // Verify
    double dtMax = 0;
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                auto x = x0 + ix*dx;
                auto y = y0 + iy*dy;
                auto z = z0 + iz*dz;
                auto delX = x - xSrc;
                auto delY = y - ySrc;
                auto delZ = z - zSrc;
                auto t = std::sqrt(delX*delX + delY*delY + delZ*delZ)/velocity;
                auto indx = iz*nx*ny + iy*nx + ix;
                dtMax = std::max(dtMax, std::abs(t - travelTimes.at(indx)));
            }
        }
    }
    EXPECT_NEAR(dtMax, 0, 1.e-15);
}

TEST(Analytic, LinearGradient2D)
{
    std::pair<double, double> velocity({1000, 4000});
    int nx = 56; 
    int nz = 24; 
    double dx = 200;
    double dz = 201;
    double x0 = 0;
    double z0 = 1;
    double xSrc = (nx - 1)/4.*dx + x0 + 4;
    double zSrc = (nz - 1)/4.*dz + z0 - 40;   
    // Create geometry
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
    // Initialize solver
    EikonalXX::Analytic::LinearGradient2D<double> solver;
    EXPECT_NO_THROW(solver.initialize(geometry));
    EXPECT_NO_THROW(solver.setSource(std::pair(xSrc, zSrc)));
    EXPECT_NO_THROW(solver.setVelocityModel(velocity));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_TRUE(solver.haveSource());
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_NO_THROW(solver.solve());
    auto travelTimes = solver.getTravelTimeField();
solver.writeVTK("test_lg.vtk");
}

}
