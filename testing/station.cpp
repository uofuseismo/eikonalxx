#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/station2d.hpp"
#include "eikonalxx/station3d.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Station, station2d)
{
    int nx = 21;
    int nz = 20;
    double dx = 100; 
    double dz = 100.05;
    double x0 = 100;
    double z0 = 200; 
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double zSrc = z0 + dz*(nz/2);
    std::string stationName = "UU.FORK";

    // Initialize the geometry
    Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);

    // Set the station
    Station2D station;
    station.setName(stationName);
    EXPECT_NO_THROW(station.setGeometry(geometry));
    EXPECT_NO_THROW(station.setLocationInX(xSrc));
    EXPECT_NO_THROW(station.setLocationInZ(zSrc));
    // Check copy c'tor + assignment
    Station2D stationCopy(station);
    // Check geometry was correctly set
    auto geoCopy = station.getGeometry();
    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-14);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-14);
    // And now the station parameters
    EXPECT_EQ(station.getName(), stationName);
    EXPECT_NEAR(stationCopy.getOffsetInX(), xSrc - x0, 1.e-14);
    EXPECT_NEAR(stationCopy.getOffsetInZ(), zSrc - z0, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInX(), xSrc, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInZ(), zSrc, 1.e-14);
    EXPECT_EQ(stationCopy.getCellInX(), 10);
    EXPECT_EQ(stationCopy.getCellInZ(), 10);
    EXPECT_EQ(stationCopy.getCell(), (nx - 1)*10 + 10); 
    // Shift station to free surface
    stationCopy.setZToFreeSurface(); 
    EXPECT_NEAR(stationCopy.getOffsetInZ(), 0, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInZ(), z0, 1.e-14);
    EXPECT_EQ(stationCopy.getCellInZ(), 0);
}

TEST(Station, station3d)
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
    std::string stationName = "UU.HRU";

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

    // Set the station
    Station3D station;
    station.setName(stationName);
    EXPECT_NO_THROW(station.setGeometry(geometry));
    EXPECT_NO_THROW(station.setLocationInX(xSrc));
    EXPECT_NO_THROW(station.setLocationInY(ySrc));
    EXPECT_NO_THROW(station.setLocationInZ(zSrc));

    // Check copy c'tor + assignment
    Station3D stationCopy(station);
    // Check geometry was correctly set
    auto geoCopy = station.getGeometry();
    EXPECT_NEAR(geoCopy.getGridSpacingInX(), dx, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInY(), dy, 1.e-14);
    EXPECT_NEAR(geoCopy.getGridSpacingInZ(), dz, 1.e-14);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInX(), nx);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInY(), ny);
    EXPECT_EQ(geoCopy.getNumberOfGridPointsInZ(), nz);
    EXPECT_NEAR(geoCopy.getOriginInX(), x0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInY(), y0, 1.e-14);
    EXPECT_NEAR(geoCopy.getOriginInZ(), z0, 1.e-14);
    // And now the station parameters
    EXPECT_EQ(station.getName(), stationName);
    EXPECT_NEAR(stationCopy.getOffsetInX(), xSrc - x0, 1.e-14);
    EXPECT_NEAR(stationCopy.getOffsetInY(), ySrc - y0, 1.e-14);
    EXPECT_NEAR(stationCopy.getOffsetInZ(), zSrc - z0, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInX(), xSrc, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInY(), ySrc, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInZ(), zSrc, 1.e-14);
    EXPECT_EQ(stationCopy.getCellInX(), 10);
    EXPECT_EQ(stationCopy.getCellInY(),  8);
    EXPECT_EQ(stationCopy.getCellInZ(), 10);
    EXPECT_EQ(stationCopy.getCell(),
              10*(nx - 1)*(ny - 1) + 8*(nx - 1) + 10);
    // Shift station to free surface
    stationCopy.setZToFreeSurface(); 
    EXPECT_NEAR(stationCopy.getOffsetInZ(), 0, 1.e-14);
    EXPECT_NEAR(stationCopy.getLocationInZ(), z0, 1.e-14);
    EXPECT_EQ(stationCopy.getCellInZ(), 0); 
}


}
