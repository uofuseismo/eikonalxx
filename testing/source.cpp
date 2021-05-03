#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/source3d.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Source, source2d)
{
    int nx = 21;
    int nz = 20;
    double dx = 100; 
    double dz = 100.05;
    double x0 = 100;
    double z0 = 200; 
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double zSrc = z0 + dz*(nz/2);

    // Initialize the geometry
    Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);

    // Set the source
    Source2D source;
    EXPECT_NO_THROW(source.setGeometry(geometry));
    EXPECT_NO_THROW(source.setLocationInX(xSrc));
    EXPECT_NO_THROW(source.setLocationInZ(zSrc));
    // Check copy c'tor + assignment
    Source2D sourceCopy(source);
    // Check geometry was correctly set
    auto geoCopy = source.getGeometry();
    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-14);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-14);
    // And now the source parameters
    EXPECT_NEAR(sourceCopy.getOffsetInX(), xSrc - x0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getOffsetInZ(), zSrc - z0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInX(), xSrc, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInZ(), zSrc, 1.e-14);
    EXPECT_EQ(sourceCopy.getCellInX(), 10);
    EXPECT_EQ(sourceCopy.getCellInZ(), 10);
    EXPECT_EQ(sourceCopy.getCell(), (nx - 1)*10 + 10); 
    // Shift source to free surface
    sourceCopy.setZToFreeSurface(); 
    EXPECT_NEAR(sourceCopy.getOffsetInZ(), 0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInZ(), z0, 1.e-14);
    EXPECT_EQ(sourceCopy.getCellInZ(), 0);
}

TEST(Source, source3d)
{
    int nx = 21; 
    int ny = 18;
    int nz = 20; 
    double dx = 100; 
    double dy = 200;
    double dz = 100.05;
    double x0 = 100;
    double y0 =-10;
    double z0 = 200; 
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double ySrc = y0 + dy*(ny/2) - dy/4;
    double zSrc = z0 + dz*(nz/2);

    // Initialize the geometry
    Geometry3D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInY(ny);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInY(dy);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInY(y0);
    geometry.setOriginInZ(z0);

    // Set the source
    Source3D source;
    EXPECT_NO_THROW(source.setGeometry(geometry));
    EXPECT_NO_THROW(source.setLocationInX(xSrc));
    EXPECT_NO_THROW(source.setLocationInY(ySrc));
    EXPECT_NO_THROW(source.setLocationInZ(zSrc));

    // Check copy c'tor + assignment
    Source3D sourceCopy(source);
    // Check geometry was correctly set
    auto geoCopy = source.getGeometry();
    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInY(), dy, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-14);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInY(), ny);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInY(), y0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-14);
    // And now the source parameters
    EXPECT_NEAR(sourceCopy.getOffsetInX(), xSrc - x0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getOffsetInY(), ySrc - y0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getOffsetInZ(), zSrc - z0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInX(), xSrc, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInY(), ySrc, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInZ(), zSrc, 1.e-14);
    EXPECT_EQ(sourceCopy.getCellInX(), 10);
    EXPECT_EQ(sourceCopy.getCellInY(),  8);
    EXPECT_EQ(sourceCopy.getCellInZ(), 10);
    EXPECT_EQ(sourceCopy.getCell(),
              10*(nx - 1)*(ny - 1) + 8*(nx - 1) + 10);
    // Shift source to free surface
    sourceCopy.setZToFreeSurface(); 
    EXPECT_NEAR(sourceCopy.getOffsetInZ(), 0, 1.e-14);
    EXPECT_NEAR(sourceCopy.getLocationInZ(), z0, 1.e-14);
    EXPECT_EQ(sourceCopy.getCellInZ(), 0); 
}


}
