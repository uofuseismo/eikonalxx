#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/ray/gradientTracer2d.hpp"
#include "eikonalxx/ray/gradientTracerOptions.hpp"
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/path3d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/point3d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/segment3d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/source3d.hpp"
#include "eikonalxx/station2d.hpp"
#include "eikonalxx/station3d.hpp"
#include "eikonalxx/analytic/homogeneous2d.hpp"
#include "eikonalxx/analytic/homogeneous3d.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Ray, Point2D)
{
    const double x{5};
    const double z{7};

    Ray::Point2D point{x, z};
    //point.setPositionInX(x);
    //point.setPositionInZ(z);

    Ray::Point2D copy(point);
    EXPECT_TRUE(copy.havePositionInX());
    EXPECT_TRUE(copy.havePositionInZ());
    EXPECT_NEAR(copy.getPositionInX(), x, 1.e-14);
    EXPECT_NEAR(copy.getPositionInZ(), z, 1.e-14);

    point.clear();
    EXPECT_FALSE(point.havePositionInX());
    EXPECT_FALSE(point.havePositionInZ());
}

TEST(Ray, Point3D)
{
    const double x{5};
    const double y{6};
    const double z{7};

    Ray::Point3D point{x, y, z};
    //point.setPositionInX(x);
    //point.setPositionInY(y);
    //point.setPositionInZ(z);

    Ray::Point3D copy(point);
    EXPECT_TRUE(copy.havePositionInX());
    EXPECT_TRUE(copy.havePositionInY());
    EXPECT_TRUE(copy.havePositionInZ());
    EXPECT_NEAR(copy.getPositionInX(), x, 1.e-14);
    EXPECT_NEAR(copy.getPositionInY(), y, 1.e-14);
    EXPECT_NEAR(copy.getPositionInZ(), z, 1.e-14);

    point.clear();
    EXPECT_FALSE(point.havePositionInX());
    EXPECT_FALSE(point.havePositionInY());
    EXPECT_FALSE(point.havePositionInZ());
}

TEST(Ray, Segment2D)
{
    const double velocity{3500};
    const double x0{100};
    const double z0{200};
    const double x1{450};
    const double z1{675};
    const double length = std::sqrt( (x1 - x0)*(x1 - x0)
                                   + (z1 - z0)*(z1 - z0) );
    const double travelTime = length/velocity;
    int cellIndex = 254;

    Ray::Point2D start;
    start.setPositionInX(x0);
    start.setPositionInZ(z0);

    Ray::Point2D end;
    end.setPositionInX(x1);
    end.setPositionInZ(z1);

    Ray::Segment2D segment;
    EXPECT_NO_THROW(segment.setStartAndEndPoint(std::pair {start, end}));
    EXPECT_NO_THROW(segment.setVelocity(velocity));
    EXPECT_NO_THROW(segment.setVelocityModelCellIndex(cellIndex));

    Ray::Segment2D copy(segment);
    EXPECT_NEAR(copy.getStartPoint().getPositionInX(), x0, 1.e-15);
    EXPECT_NEAR(copy.getStartPoint().getPositionInZ(), z0, 1.e-15);
    EXPECT_NEAR(copy.getEndPoint().getPositionInX(), x1, 1.e-15);
    EXPECT_NEAR(copy.getEndPoint().getPositionInZ(), z1, 1.e-15);
    EXPECT_NEAR(copy.getVelocity(), velocity, 1.e-12);
    EXPECT_NEAR(copy.getSlowness(), 1./velocity, 1.e-15);
    EXPECT_NEAR(copy.getLength(),   length, 1.e-15);
    EXPECT_NEAR(copy.getTravelTime(), travelTime, 1.e-15);
    EXPECT_EQ(copy.getVelocityModelCellIndex(), cellIndex);

    segment.clear();
    EXPECT_FALSE(segment.haveStartAndEndPoint());
    EXPECT_FALSE(segment.haveVelocity());
    EXPECT_FALSE(segment.haveVelocityModelCellIndex());
}

TEST(Ray, Segment3D)
{
    const double velocity{3500};
    const double x0{100};
    const double y0{150};
    const double z0{200};
    const double x1{450};
    const double y1{505};
    const double z1{675};
    const double length = std::sqrt( (x1 - x0)*(x1 - x0) 
                                   + (y1 - y0)*(y1 - y0)
                                   + (z1 - z0)*(z1 - z0) );
    const double travelTime = length/velocity;
    int cellIndex = 634;

    Ray::Point3D start;
    start.setPositionInX(x0);
    start.setPositionInY(y0);
    start.setPositionInZ(z0);

    Ray::Point3D end;
    end.setPositionInX(x1);
    end.setPositionInY(y1);
    end.setPositionInZ(z1);

    Ray::Segment3D segment;
    EXPECT_NO_THROW(segment.setStartAndEndPoint(std::pair {start, end}));
    EXPECT_NO_THROW(segment.setVelocity(velocity));
    EXPECT_NO_THROW(segment.setVelocityModelCellIndex(cellIndex));

    Ray::Segment3D copy(segment);
    EXPECT_NEAR(copy.getStartPoint().getPositionInX(), x0, 1.e-15);
    EXPECT_NEAR(copy.getStartPoint().getPositionInY(), y0, 1.e-15);
    EXPECT_NEAR(copy.getStartPoint().getPositionInZ(), z0, 1.e-15);
    EXPECT_NEAR(copy.getEndPoint().getPositionInX(), x1, 1.e-15);
    EXPECT_NEAR(copy.getEndPoint().getPositionInY(), y1, 1.e-15);
    EXPECT_NEAR(copy.getEndPoint().getPositionInZ(), z1, 1.e-15);
    EXPECT_NEAR(copy.getVelocity(), velocity, 1.e-12);
    EXPECT_NEAR(copy.getSlowness(), 1./velocity, 1.e-15);
    EXPECT_NEAR(copy.getLength(),   length, 1.e-15);
    EXPECT_NEAR(copy.getTravelTime(), travelTime, 1.e-15);
    EXPECT_EQ(copy.getVelocityModelCellIndex(), cellIndex);

    segment.clear();
    EXPECT_FALSE(segment.haveStartAndEndPoint());
    EXPECT_FALSE(segment.haveVelocity());
    EXPECT_FALSE(segment.haveVelocityModelCellIndex());
}

TEST(Ray, Path2D)
{
    const std::vector<double> xs{0,   100, 200, 300, 400};
    const std::vector<double> zs{150, 250, 300, 320, 333};
    const std::vector<double> vels{500, 600, 700, 800};
    std::vector<Ray::Segment2D> segments;
    double travelTime = 0;
    double length = 0;
    for (int i = 0; i < static_cast<int> (xs.size()) - 1; ++i)
    {
        Ray::Point2D point0, point1;
        auto dx = xs[i + 1] - xs[i];
        auto dz = zs[i + 1] - zs[i];
        travelTime = travelTime + std::sqrt(dx*dx + dz*dz)/vels[i];
        length = length + std::sqrt(dx*dx + dz*dz);
        point0.setPositionInX(xs[i]);
        point0.setPositionInZ(zs[i]); 
        point1.setPositionInX(xs[i + 1]);
        point1.setPositionInZ(zs[i + 1]);
        Ray::Segment2D segment;
        segment.setStartAndEndPoint(std::pair {point0, point1});
        EXPECT_NO_THROW(segment.setVelocity(vels.at(i)));
        segments.push_back(std::move(segment));
    }
 
    Ray::Path2D path; 
    path.open();
    EXPECT_TRUE(path.isOpen());
    for (const auto &s : segments)
    {
        EXPECT_NO_THROW(path.append(s));
    }
    EXPECT_NO_THROW(path.close());
    EXPECT_FALSE(path.isOpen());
    EXPECT_EQ(path.size(), segments.size()); 
    EXPECT_NEAR(path.getTravelTime(), travelTime, 1.e-10);
    EXPECT_NEAR(path.getLength(), length, 1.e-7);

    EXPECT_NO_THROW(path.set(segments));
    path.clear();
    EXPECT_FALSE(path.isOpen());
    EXPECT_NO_THROW(path.set(segments));
    EXPECT_FALSE(path.isOpen());

    Ray::Path2D copy(path);
    EXPECT_NEAR(copy.getTravelTime(), travelTime, 1.e-10);
    EXPECT_NEAR(copy.getLength(), length, 1.e-7);
    int is = 0;
    for (const auto &s : copy)
    {
        EXPECT_NEAR(s.getStartPoint().getPositionInX(), 
                    segments[is].getStartPoint().getPositionInX(), 1.e-15);
        EXPECT_NEAR(s.getStartPoint().getPositionInZ(),
                    segments[is].getStartPoint().getPositionInZ(), 1.e-15);
        EXPECT_NEAR(s.getEndPoint().getPositionInX(),
                    segments[is].getEndPoint().getPositionInX(), 1.e-15);
        EXPECT_NEAR(s.getEndPoint().getPositionInZ(), 
                    segments[is].getEndPoint().getPositionInZ(), 1.e-15);
        EXPECT_NEAR(s.getVelocity(),
                    segments[is].getVelocity(), 1.e-15);
        is = is + 1;
    }
}

TEST(Ray, Path3D)
{
    const std::vector<double> xs{0,   100, 200, 300, 400};
    const std::vector<double> ys{24,   32,  84,  96, 101};
    const std::vector<double> zs{150, 250, 300, 320, 333};
    const std::vector<double> vels{500, 600, 700, 800};
    std::vector<Ray::Segment3D> segments;
    double travelTime = 0;
    double length = 0;
    for (int i = 0; i < static_cast<int> (xs.size()) - 1; ++i)
    {   
        Ray::Point3D point0, point1;
        auto dx = xs[i + 1] - xs[i];
        auto dy = ys[i + 1] - ys[i];
        auto dz = zs[i + 1] - zs[i];
        travelTime = travelTime + std::sqrt(dx*dx + dy*dy + dz*dz)/vels[i];
        length = length + std::sqrt(dx*dx + dy*dy + dz*dz);
        point0.setPositionInX(xs[i]);
        point0.setPositionInY(ys[i]); 
        point0.setPositionInZ(zs[i]); 
        point1.setPositionInX(xs[i + 1]);
        point1.setPositionInY(ys[i + 1]);
        point1.setPositionInZ(zs[i + 1]);
        Ray::Segment3D segment;
        segment.setStartAndEndPoint(std::pair {point0, point1});
        EXPECT_NO_THROW(segment.setVelocity(vels.at(i)));
        segments.push_back(std::move(segment));
    } 

    Ray::Path3D path;
    path.open();
    EXPECT_TRUE(path.isOpen());
    for (const auto &s : segments)
    {
        EXPECT_NO_THROW(path.append(s));
    }
    EXPECT_NO_THROW(path.close());
    EXPECT_FALSE(path.isOpen());
    EXPECT_EQ(path.size(), segments.size());
    EXPECT_NEAR(path.getTravelTime(), travelTime, 1.e-10);
    EXPECT_NEAR(path.getLength(), length, 1.e-7);

    EXPECT_NO_THROW(path.set(segments));
    path.clear();
    EXPECT_FALSE(path.isOpen());
    EXPECT_NO_THROW(path.set(segments));
    EXPECT_FALSE(path.isOpen());

    Ray::Path3D copy(path);
    EXPECT_NEAR(copy.getTravelTime(), travelTime, 1.e-10);
    EXPECT_NEAR(copy.getLength(), length, 1.e-7);
    int is = 0;
    for (const auto &s : copy)
    {
        EXPECT_NEAR(s.getStartPoint().getPositionInX(),
                    segments[is].getStartPoint().getPositionInX(), 1.e-15);
        EXPECT_NEAR(s.getStartPoint().getPositionInY(),
                    segments[is].getStartPoint().getPositionInY(), 1.e-15);
        EXPECT_NEAR(s.getStartPoint().getPositionInZ(),
                    segments[is].getStartPoint().getPositionInZ(), 1.e-15);
        EXPECT_NEAR(s.getEndPoint().getPositionInX(),
                    segments[is].getEndPoint().getPositionInX(), 1.e-15);
        EXPECT_NEAR(s.getEndPoint().getPositionInY(),
                    segments[is].getEndPoint().getPositionInY(), 1.e-15);
        EXPECT_NEAR(s.getEndPoint().getPositionInZ(),
                    segments[is].getEndPoint().getPositionInZ(), 1.e-15);
        EXPECT_NEAR(s.getVelocity(),
                    segments[is].getVelocity(), 1.e-15);
        is = is + 1;
    }
}

TEST(Ray, GradientTracerOptions)
{
    EikonalXX::Ray::GradientTracerOptions options;
    std::vector<std::pair<int, double>> radiusScaleFactor{ std::pair {2, 0.2},
                                                           std::pair {1, 0.1} };
    std::vector<std::pair<int, double>> reference{ std::pair {1, 0.1},
                                                   std::pair {2, 0.2},
                                                   std::pair {std::numeric_limits<int>::max(), 0.2} };
    options.setRadiusScaleFactor(radiusScaleFactor);
 
    EikonalXX::Ray::GradientTracerOptions copy(options);
    auto rsf = copy.getRadiusScaleFactor();
    EXPECT_EQ(rsf.size(), reference.size());
    for (int i = 0; i < static_cast<int> (rsf.size()); ++i)
    {
        EXPECT_EQ(rsf[i].first, reference[i].first);
        EXPECT_NEAR(rsf[i].second, reference[i].second, 1.e-14);
    }
}

TEST(Ray, HomogeneousGradientTracer2D)
{
    constexpr double velocity{4000};
    const double dx{100};
    const double dz{101};
    const double x0{1};
    const double z0{2};
    const int nx{201};
    const int nz{101};
    const double xs{x0 + nx/2*dx};
    const double zs{z0 + nz/2*dz};
    const int nStationsVertical{11};
    const int nStationsHorizontal{11};
    // Create geometry
    EikonalXX::Geometry2D geometry;
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
    // Create source
    EikonalXX::Source2D source;
    source.setGeometry(geometry);
    source.setLocationInX(xs);
    source.setLocationInZ(zs);
    // Get the stations
    std::vector<EikonalXX::Station2D> stations;
    EikonalXX::Station2D station;
    station.setGeometry(geometry);
    station.setName("UU.T" + std::to_string(stations.size()));
    station.setLocationInX(xs + 0.001);
    station.setLocationInZ(zs + 0.001);
    stations.push_back(station);
    // Distribute set of stations along free surface
    for (int is = 0; is < nStationsHorizontal; ++is)
    {
        station.setName("UU.T" + std::to_string(stations.size()));
        station.setLocationInX(x0 + is/(nStationsHorizontal - 1.0)*(nx - 1)*dx);
        station.setZToFreeSurface();
        stations.push_back(station);
    }
    // Distribute set of stations through depth
    for (int is = 0; is < nStationsVertical; ++is)
    {
        station.setName("UU.T" + std::to_string(stations.size()));
        station.setLocationInX(x0);
        station.setLocationInZ(z0 + is/(nStationsVertical - 1.0)*(nz - 1)*dz);
        stations.push_back(station);

        station.setName("UU.T" + std::to_string(stations.size()));
        station.setLocationInX(x0 + (nx - 1)*dx);
        stations.push_back(station);
    }
    // Solve eikonal equation and compute gradient
    EikonalXX::Analytic::Homogeneous2D<double> solver;
    solver.initialize(geometry);
    solver.setVelocityModel(velocity);
    solver.setSource(source);
    solver.solve();
    solver.computeTravelTimeGradientField();
    // Create the ray tracer
    EikonalXX::Ray::GradientTracerOptions options;
    EikonalXX::Ray::GradientTracer2D tracer;
    EXPECT_NO_THROW(tracer.initialize(options, geometry));
    EXPECT_NO_THROW(tracer.setStations(stations));
    tracer.trace(solver);
solver.writeVTK("homog2d.vtk");
tracer.writeVTK("raypaths.vtk");
}

}
