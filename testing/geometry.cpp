#include <iostream>
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(TestGeometry, geometry2d)
{
    Geometry2D geo;
    int nx = 5;
    int nz = 6;
    double dx = 3;
    double dz = 4;
    double x0 = 1;
    double z0 = 2;

    EXPECT_NO_THROW(geo.setNumberOfGridPointsInX(nx));
    EXPECT_NO_THROW(geo.setNumberOfGridPointsInZ(nz));
    EXPECT_NO_THROW(geo.setGridSpacingInX(dx));
    EXPECT_NO_THROW(geo.setGridSpacingInZ(dz));
    geo.setOriginInX(x0);
    geo.setOriginInZ(z0);

    // Use copy c'tor and check equality
    Geometry2D geoCopy(geo);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_EQ(geoCopy.getNumberOfCellsInX(), nx - 1);
    EXPECT_EQ(geoCopy.getNumberOfCellsInZ(), nz - 1);
    EXPECT_EQ(geoCopy.getNumberOfGridPoints(), nx*nz);
    EXPECT_EQ(geoCopy.getNumberOfCells(), (nx - 1)*(nz - 1));

    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-10);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-10);

    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-10);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-10);
}

TEST(TestGeometry, geometry3d)
{
    Geometry3D geo;
    int nx = 5;
    int ny = 7;
    int nz = 6;
    double dx = 3;
    double dy = 2;
    double dz = 4;
    double x0 = 1;
    double y0 =-1;
    double z0 = 2;

    EXPECT_NO_THROW(geo.setNumberOfGridPointsInX(nx));
    EXPECT_NO_THROW(geo.setNumberOfGridPointsInY(ny));
    EXPECT_NO_THROW(geo.setNumberOfGridPointsInZ(nz));
    EXPECT_NO_THROW(geo.setGridSpacingInX(dx));
    EXPECT_NO_THROW(geo.setGridSpacingInY(dy));
    EXPECT_NO_THROW(geo.setGridSpacingInZ(dz));
    geo.setOriginInX(x0);
    geo.setOriginInY(y0);
    geo.setOriginInZ(z0);

    // Use copy c'tor and check equality
    Geometry3D geoCopy(geo);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInY(), ny);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_EQ(geoCopy.getNumberOfCellsInX(), nx - 1); 
    EXPECT_EQ(geoCopy.getNumberOfCellsInY(), ny - 1);
    EXPECT_EQ(geoCopy.getNumberOfCellsInZ(), nz - 1); 
    EXPECT_EQ(geoCopy.getNumberOfGridPoints(), nx*ny*nz);
    EXPECT_EQ(geoCopy.getNumberOfCells(), (nx - 1)*(ny - 1)*(nz - 1));

    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-10);
    EXPECT_NEAR(geoCopy.getGridSpacingInY(), dy, 1.e-10);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-10);

    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-10);
    EXPECT_NEAR(geoCopy.getOriginInY(), y0, 1.e-10);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-10);
}

}
