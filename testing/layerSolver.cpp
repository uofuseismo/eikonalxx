#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/layerSolver.hpp"
#include "eikonalxx/ray/gradientTracer2d.hpp"
#include "eikonalxx/ray/gradientTracerOptions.hpp"
#include "eikonalxx/solver2d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/station2d.hpp"
#include "ray/layerTracer.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX::Ray;

bool operator==(const Point2D &lhs, const Point2D &rhs)
{
    if (std::abs(lhs.getPositionInX() - rhs.getPositionInX()) > 1.e-1)
    {
        return false;
    }
    if (std::abs(lhs.getPositionInZ() - rhs.getPositionInZ()) > 1.e-1)
    {
        return false;
    }
    return true;
}

bool operator==(const Segment2D &lhs, const Segment2D &rhs)
{
    if (lhs.getVelocityModelCellIndex() != rhs.getVelocityModelCellIndex())
    {
        std::cerr << "Segment cell model index not equal" << std::endl;
        return false;
    }
    if (std::abs(lhs.getLength() - rhs.getLength()) > 1.e-1)
    {
        std::cerr << "Segment length not equal" << std::endl;
        return false;
    }
    if (std::abs(lhs.getTravelTime() - rhs.getTravelTime()) > 1.e-1)
    {
        std::cerr << "Segment travel time not equal" << std::endl;
        return false;
    }
    if (lhs.getStartPoint() != rhs.getStartPoint())
    {
        std::cerr << "Segment start point not equal" << std::endl;
        return false;
    }
    if (lhs.getEndPoint() != rhs.getEndPoint())
    {
        std::cerr << "Segment end point not equal" << std::endl;
        return false;
    }
    return true;
}

bool operator==(const Path2D &lhs, const Path2D &rhs)
{
    if (lhs.size() != rhs.size())
    {
        std::cerr << "Number of segments differs" << std::endl;
        return false;
    }
    if (std::abs(lhs.getTravelTime() - rhs.getTravelTime()) > 1.e-6)
    {
        std::cerr << "Path travel time differs: "
                  << lhs.getTravelTime() << "," << rhs.getTravelTime()
                  << std::endl;
        return false;
    }
    if (std::abs(lhs.getLength() - rhs.getLength()) > 1.e-1)
    {
        std::cerr << "Path lengths differ: "
                  << lhs.getLength() << "," << rhs.getLength() << std::endl;
        return false;
    }
    for (int i = 0; i < static_cast<int> (lhs.size()); ++i)
    {
        if (lhs.at(i) != rhs.at(i)){return false;}
    }
    return true;
}

bool operator==(const std::vector<Path2D> &lhs, const std::vector<Path2D> &rhs)
{
    if (lhs.size() != rhs.size())
    {
        std::cerr << "Number of ray paths differs" << std::endl;
        return false;
    }
    for (size_t i = 0; i < lhs.size(); ++i)
    {
        if (lhs[i] != rhs[i]){return false;}
    }
    return true;
}

[[maybe_unused]]
void createReferenceSolution(const double sourceDepth,
                             const double stationDepth,
                             const std::vector<double> &interfaces,
                             const std::vector<double> &velocities,
                             const std::vector<double> &stationLocations,
                             const double dx = 50,
                             const double xWidth = 200000,
                             const double zWidth =  60000)
{
    auto dz = dx;
    auto z0 = interfaces[0];
    auto x0 =-2*dx;
    auto x1 = xWidth + 2*dx;
    auto nx = static_cast<int> ( (x1 - x0)/dx ) + 1;
    auto nz = static_cast<int> ( zWidth/dx ) + 1;
    EikonalXX::Geometry2D geometry;
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dx);
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);
  
    std::vector<double> velocityModel(nx*nz);
    for (int iz = 0; iz < nz; ++iz)
    {
        auto depth = z0 + iz*dz; 
        auto layer = std::distance(interfaces.begin(),
                                   std::upper_bound(interfaces.begin(),
                                                    interfaces.end(),
                                                    depth)) - 1;
        auto velocity = velocities.at(layer);
        auto i0 = iz*nx;
        auto i1 = (iz + 1)*nx;
        std::fill(velocityModel.data() + i0,
                  velocityModel.data() + i1,
                  velocity);
    }
    std::cout << "Setting velocity model..." << std::endl;
    EikonalXX::Model2D<double> model;
    model.initialize(geometry);
    model.setNodalVelocities(velocityModel.size(),
                             velocityModel.data(),
                             EikonalXX::Ordering2D::Natural);

    EikonalXX::Source2D source;
    source.setGeometry(geometry);
    source.setLocationInX(0);
    source.setLocationInZ(sourceDepth);

    std::vector<EikonalXX::Station2D> stations;
    for (const auto &stationLocation : stationLocations)
    {
        EikonalXX::Station2D station;
        station.setGeometry(geometry);
        station.setLocationInX(stationLocation);
        station.setLocationInZ(stationDepth);
        stations.push_back(station);
    }  

    std::cout << "Initializing solver..." << std::endl;
    EikonalXX::Solver2D<double> solver;
    EikonalXX::SolverOptions options;
    options.setFactoredEikonalEquationSolverRadius(15);
    options.setNumberOfSweeps(4);
    options.setConvergenceTolerance(1.e-14);
    options.setAlgorithm(EikonalXX::SolverAlgorithm::FastSweepingMethod);

    solver.initialize(options, geometry);
    solver.setVelocityModel(model);
    solver.setSource(source);
    solver.setStations(stations);
    std::cout << "Solving..." << std::endl;
    solver.solve(); 
    std::cout << "Compute gradient..." << std::endl;
    solver.computeTravelTimeGradientField();
    solver.writeVTK("uussTTimes.vtk", "travel_time_field_s", true);

    auto travelTimes = solver.getTravelTimesAtStations();
    for (const auto &travelTime : travelTimes)
    {
        std::cout << std::setprecision(12) << "Eikonal travel times: " << travelTime << std::endl;
    }

    EikonalXX::Ray::GradientTracerOptions tracerOptions;
    EikonalXX::Ray::GradientTracer2D gradTracer;
    gradTracer.initialize(tracerOptions, geometry);
    gradTracer.setStations(stations);
    gradTracer.trace(solver);
    gradTracer.writeVTK("uussRayPaths.vtk");
    auto rayPaths = gradTracer.getRayPaths();
    for (const auto &rayPath : rayPaths)
    {
        std::cout << "Take off angle: " << rayPath.getTakeOffAngle() << "," << rayPath.getTravelTime() << std::endl;
    }
}

TEST(Ray, LayerSolverHalfSpaceInternal)
{
    const double velocity{3500};
    const double sourceDepth{0};
    const double stationOffset{5000};
    const double stationDepth{-4500};
    EikonalXX::Ray::Path2D path;
    auto returnCode = ::traceWholeSpace(1./velocity,
                                        sourceDepth,
                                        stationOffset,
                                        stationDepth,
                                        &path);
    EXPECT_EQ(returnCode, ::ReturnCode::Hit); 

    Point2D point0{0, 0}; 
    Point2D point1{5000, -4500};
    Segment2D segment;
    segment.setStartAndEndPoint(std::pair {point0, point1}); 
    segment.setVelocity(velocity);
    segment.setVelocityModelCellIndex(0);
    Path2D referencePath;
    referencePath.open();
    referencePath.append(segment);
    referencePath.close();
    EXPECT_TRUE(referencePath == path);
}

TEST(Ray, LayerSolverVerticalReflectionInternal)
{
    constexpr int nLayers = 5;
    constexpr double convergence{1}; // meters
    auto interfaces = ::augmentInterfacesVector({-4500,   40,  15600, 26500, 40500});
    auto velocities = ::augmentVelocityVector({       3500, 5900,  6400,  7500,  7900});
    auto slownesses = ::toSlownessVector(velocities);
    const double sourceDepth{20000};
    const double stationDepth{-4000};
    const double stationOffset{0};
    constexpr auto isAugmented{true};
    auto sourceLayer  = ::getLayer(sourceDepth,  interfaces, isAugmented);
    auto stationLayer = ::getLayer(stationDepth, interfaces, isAugmented);
    EXPECT_EQ(sourceLayer,  2);
    EXPECT_EQ(stationLayer, 0); 
    // Draw the reference rays
    Point2D sourcePoint{0, sourceDepth};
    Point2D stationPoint{0, stationDepth};
    Point2D interface1{0, 40};
    Point2D interface2{0, 15600};
    Point2D interface3{0, 26500};
    Point2D interface4{0, 40500};  

    Segment2D segmentDown23;
    segmentDown23.setStartAndEndPoint(std::pair {sourcePoint, interface3});
    segmentDown23.setVelocity(6400);
    segmentDown23.setVelocityModelCellIndex(2);

    Segment2D segmentDown34;
    segmentDown34.setStartAndEndPoint(std::pair {interface3, interface4});
    segmentDown34.setVelocity(7500);
    segmentDown34.setVelocityModelCellIndex(3);

    auto segmentUp43 = segmentDown34;
    segmentUp43.reverse();

    Segment2D segmentUp32;
    segmentUp32.setStartAndEndPoint(std::pair {interface3, interface2});
    segmentUp32.setVelocity(6400);
    segmentUp32.setVelocityModelCellIndex(2);

    Segment2D segmentUp21;
    segmentUp21.setStartAndEndPoint(std::pair {interface2, interface1});
    segmentUp21.setVelocity(5900);
    segmentUp21.setVelocityModelCellIndex(1);

    Segment2D segmentUp10;
    segmentUp10.setStartAndEndPoint(std::pair {interface1, stationPoint});
    segmentUp10.setVelocity(3500);
    segmentUp10.setVelocityModelCellIndex(0);

    std::vector<Segment2D> rayPath1{segmentDown23,
                                    segmentUp32, segmentUp21, segmentUp10};
    std::vector<Segment2D> rayPath2{segmentDown23, segmentDown34,
                                    segmentUp43, segmentUp32, segmentUp21, segmentUp10};
    Path2D path1;
    path1.set(rayPath1);
    Path2D path2;
    path2.set(rayPath2);
    std::vector<Path2D> rayPaths{path1, path2};
    // Shoot down
    int iPath = 0;
    for (int layer = sourceLayer; layer < nLayers - 1; ++layer)
    {
        std::vector<::Segment> segments;
        auto returnCode = ::traceVerticalReflectionDown(interfaces,
                                      slownesses,
                                      sourceLayer,
                                      layer,
                                      sourceDepth,
                                      stationLayer,
                                      stationDepth,
                                      stationOffset,
                                      &segments,
                                      convergence);
        EXPECT_EQ(returnCode, ReturnCode::Hit);
        auto rayPath = ::toRayPath(segments);
        EXPECT_TRUE(rayPath == rayPaths[iPath]);
        iPath = iPath + 1;
    }

    // Shoot down to station
    iPath = 0;
    for (int layer = sourceLayer; layer < nLayers - 1; ++layer) 
    {
        std::vector<::Segment> segments;
        auto returnCode = ::traceVerticalReflectionDown(interfaces,
                                                        slownesses,
                                                        stationLayer,
                                                        layer,
                                                        stationDepth,
                                                        sourceLayer,
                                                        sourceDepth,
                                                        stationOffset,
                                                        &segments,
                                                        convergence);
        EXPECT_EQ(returnCode, ReturnCode::Hit);
        auto rayPath = ::toRayPath(segments);
        auto reversedPath = rayPaths[iPath];
        reversedPath.reverse(); 
        EXPECT_TRUE(rayPath == reversedPath);
        iPath = iPath + 1;
    } 
}

TEST(Ray, LayerSolverDownGoingSameSourceStationDepthInternal)
{
    constexpr int nLayers{5};
    constexpr double depth{-2000};
    auto interfaces = ::augmentInterfacesVector({-4500,   50,  15600, 26500, 40500});
    auto velocities = ::augmentVelocityVector({       3500, 5900,  6400,  7500,  7900});
    auto slownesses = ::toSlownessVector(velocities);
    constexpr double sourceDepth{depth};
    constexpr double stationDepth{depth};
    constexpr double stationOffset{20000};
    constexpr auto isAugmented{true};
    constexpr bool allowCriticalRefractions{true};
    constexpr double hitTolerance{1};
    auto sourceLayer  = ::getLayer(sourceDepth,  interfaces, isAugmented);
    auto stationLayer = ::getLayer(stationDepth, interfaces, isAugmented);
    EXPECT_EQ(sourceLayer,  0);
    EXPECT_EQ(stationLayer, 0);

    double takeOffAngle{39.2777};

    // Shoot to first layer
    std::vector<::Segment> segments;
    auto returnCode = ::traceDownThenUp(interfaces,
                                        slownesses,
                                        takeOffAngle,
                                        sourceLayer,
                                        sourceDepth,
                                        stationLayer,
                                        stationOffset,
                                        stationDepth,
                                        0,
                                        &segments,
                                        allowCriticalRefractions,
                                        hitTolerance);
    Point2D sourcePoint{0, sourceDepth};
    Point2D point01{1676.5317885936365, 50};
    Point2D point02{18323.46821140636,  50};
    Point2D stationPoint{stationOffset, stationDepth};
    Segment2D segment01;
    Segment2D segment02;
    Segment2D segment03; 
    segment01.setStartAndEndPoint(std::pair {sourcePoint, point01});
    segment01.setVelocity(3500);
    segment01.setVelocityModelCellIndex(0);

    segment02.setStartAndEndPoint(std::pair {point01, point02});
    segment02.setVelocity(5900);
    segment02.setVelocityModelCellIndex(1);

    segment03.setStartAndEndPoint(std::pair {point02, stationPoint});
    segment03.setVelocity(3500);
    segment03.setVelocityModelCellIndex(0);

    Path2D pathTopLayer;
    pathTopLayer.set(std::vector<Segment2D> {segment01, segment02, segment03});
    EXPECT_TRUE(pathTopLayer == ::toRayPath(segments));

    // Trace at depth
    takeOffAngle = 20;
    /// Snell's law
    double theta0 = takeOffAngle*M_PI/180;
    double theta1 = std::asin((velocities[1]/velocities[0])*std::sin(theta0));
    double theta2 = std::asin((velocities[2]/velocities[1])*std::sin(theta1));
    double theta3 = std::asin((velocities[3]/velocities[2])*std::sin(theta2));
    auto dz0 = interfaces[1] - sourceDepth;//interfaces[0];
    auto dz1 = interfaces[2] - interfaces[1];
    auto dz2 = interfaces[3] - interfaces[2];
    auto dz3 = interfaces[4] - interfaces[3];

    double dx0 = dz0*std::tan(theta0); //746.1389802457148
    double dx1 = dz1*std::tan(theta1); //10972.622618039548
    double dx2 = dz2*std::tan(theta2); //8736.34515155134
    double dx3 = dz3*std::tan(theta3); //15081.670759818773
    Path2D path2;
    {
    Point2D point1{sourcePoint},
            point2{dx0, 50},
            point3{dx0 + dx0, stationDepth};
    Segment2D segment1, segment2;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[0]);
    segment1.setVelocityModelCellIndex(0); 
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[0]);
    segment2.setVelocityModelCellIndex(0); 

    path2.set( {segment1, segment2} );
    }

    Path2D path3;
    {   
    Point2D point1{sourcePoint},
            point2{dx0, 50},
            point3{dx0 + dx1, 15600.0},
            point4{dx0 + dx1 + dx1, 50},
            point5{dx0 + dx1 + dx1 + dx0, stationDepth};
    Segment2D segment1, segment2, segment3, segment4;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[0]);
    segment1.setVelocityModelCellIndex(0); 
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[1]);
    segment2.setVelocityModelCellIndex(1); 
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[1]);
    segment3.setVelocityModelCellIndex(1);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[0]);
    segment4.setVelocityModelCellIndex(0);

    path3.set( {segment1, segment2, segment3, segment4} );
    }

    Path2D path4;
    {   
    Point2D point1{sourcePoint},
            point2{dx0, 50},
            point3{dx0 + dx1, 15600.0},
            point4{dx0 + dx1 + dx2, 26500.0},
            point5{dx0 + dx1 + dx2 + dx2, 15600.0},
            point6{dx0 + dx1 + dx2 + dx2 + dx1, 50},
            point7{dx0 + dx1 + dx2 + dx2 + dx1 + dx0, stationDepth};
    Segment2D segment1, segment2, segment3, segment4, segment5, segment6;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[0]);
    segment1.setVelocityModelCellIndex(0); 
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[1]);
    segment2.setVelocityModelCellIndex(1); 
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[2]);
    segment3.setVelocityModelCellIndex(2);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[2]);
    segment4.setVelocityModelCellIndex(2);
    segment5.setStartAndEndPoint(std::pair {point5, point6});
    segment5.setSlowness(slownesses[1]);
    segment5.setVelocityModelCellIndex(1);
    segment6.setStartAndEndPoint(std::pair {point6, point7});
    segment6.setSlowness(slownesses[0]);
    segment6.setVelocityModelCellIndex(0);

    path4.set( {segment1, segment2, segment3, segment4, segment5, segment6} );
    }

    Path2D path5;
    {   
    Point2D point1{sourcePoint},
            point2{dx0, 50},
            point3{dx0 + dx1, 15600.0},
            point4{dx0 + dx1 + dx2, 26500.0},
            point5{dx0 + dx1 + dx2 + dx3, 40500.0},
            point6{dx0 + dx1 + dx2 + dx3 + dx3, 26500},
            point7{dx0 + dx1 + dx2 + dx3 + dx3 + dx2, 15600},
            point8{dx0 + dx1 + dx2 + dx3 + dx3 + dx2 + dx1, 50},
            point9{dx0 + dx1 + dx2 + dx3 + dx3 + dx2 + dx1 + dx0, stationDepth};
    Segment2D segment1, segment2, segment3, segment4, segment5, segment6, segment7, segment8;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[0]);
    segment1.setVelocityModelCellIndex(0); 
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[1]);
    segment2.setVelocityModelCellIndex(1); 
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[2]);
    segment3.setVelocityModelCellIndex(2);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[3]);
    segment4.setVelocityModelCellIndex(3);
    segment5.setStartAndEndPoint(std::pair {point5, point6});
    segment5.setSlowness(slownesses[3]);
    segment5.setVelocityModelCellIndex(3);
    segment6.setStartAndEndPoint(std::pair {point6, point7});
    segment6.setSlowness(slownesses[2]);
    segment6.setVelocityModelCellIndex(2);
    segment7.setStartAndEndPoint(std::pair {point7, point8});
    segment7.setSlowness(slownesses[1]);
    segment7.setVelocityModelCellIndex(1);
    segment8.setStartAndEndPoint(std::pair {point8, point9});
    segment8.setSlowness(slownesses[0]);
    segment8.setVelocityModelCellIndex(0);

    path5.set( {segment1, segment2, segment3, segment4, segment5, segment6, segment7, segment8} );
    }

    std::vector<Path2D> paths2{path2, path3, path4, path5};
    int iPath = 0;
    for (int layer = sourceLayer; layer < nLayers - 1; ++layer)
    {
        returnCode = ::traceDownThenUp(interfaces,
                                       slownesses,
                                       takeOffAngle,
                                       sourceLayer,
                                       sourceDepth,
                                       stationLayer,
                                       stationOffset,
                                       stationDepth,
                                       layer,
                                       &segments,
                                       true,
                                       1);
        EXPECT_TRUE(returnCode == ReturnCode::Hit ||
                    returnCode == ReturnCode::UnderShot ||
                    returnCode == ReturnCode::OverShot);
        auto rayPath = ::toRayPath(segments); 
        EXPECT_TRUE(rayPath == paths2[iPath]);
        iPath = iPath + 1;
    }
}

TEST(Ray, LayerSolverDirectInternal)
{
    //constexpr int nLayers{5};
    auto interfaces = ::augmentInterfacesVector({-4500,   50,  15600, 26500, 40500});
    auto velocities = ::augmentVelocityVector({       3500, 5900,  6400,  7500,  7900});
    auto slownesses = ::toSlownessVector(velocities);
    constexpr double sourceDepth{30000};
    constexpr double stationDepth{-2000};
    constexpr double stationOffset{9275.402998};
    constexpr auto isAugmented{true};
    auto sourceLayer  = ::getLayer(sourceDepth,  interfaces, isAugmented);
    auto stationLayer = ::getLayer(stationDepth, interfaces, isAugmented);
    EXPECT_EQ(sourceLayer,  3); 
    EXPECT_EQ(stationLayer, 0); 

    double takeOffAngle{180 - 20};
    std::vector<::Segment> segments;
    auto returnCode = ::traceDirect(interfaces,
                                    slownesses,
                                    takeOffAngle,
                                    sourceLayer,
                                    sourceDepth,
                                    stationLayer,
                                    stationOffset,
                                    stationDepth,
                                    &segments,
                                    1.0);
    EXPECT_EQ(returnCode, ReturnCode::Hit);
    auto rayPath = ::toRayPath(segments);

    auto dz3 = sourceDepth - interfaces[3];
    auto dz2 = interfaces[3] - interfaces[2];
    auto dz1 = interfaces[2] - interfaces[1];
    auto dz0 = interfaces[1] - stationDepth;

    auto theta3 = (180 - takeOffAngle)*(M_PI/180);
    auto theta2 = std::asin((velocities[2]/velocities[3])*std::sin(theta3)); 
    auto theta1 = std::asin((velocities[1]/velocities[2])*std::sin(theta2));
    auto theta0 = std::asin((velocities[0]/velocities[1])*std::sin(theta1));

    //std::cout << dz0 << "," << dz1 << "," << dz2 << "," << dz3 << std::endl;

    double dx3 = dz3*std::tan(theta3); //10972.622618039548
    double dx2 = dz2*std::tan(theta2); //8736.34515155134
    double dx1 = dz1*std::tan(theta1); //15081.670759818773
    double dx0 = dz0*std::tan(theta0);

    Path2D path;
    Point2D point1{0, sourceDepth},
            point2{dx3, interfaces[3]},
            point3{dx3 + dx2, interfaces[2]},
            point4{dx3 + dx2 + dx1, interfaces[1]},
            point5{dx3 + dx2 + dx1 + dx0, stationDepth};
    Segment2D segment1, segment2, segment3, segment4;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[3]);
    segment1.setVelocityModelCellIndex(3); 
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[2]);
    segment2.setVelocityModelCellIndex(2); 
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[1]);
    segment3.setVelocityModelCellIndex(1);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[0]);
    segment4.setVelocityModelCellIndex(0);

    path.set( {segment1, segment2, segment3, segment4} );

    EXPECT_TRUE(path == rayPath);
}

TEST(Ray, LayerSolverGeneralSourceBelowStationInternal)
{
    constexpr int nLayers{5};
    auto interfaces = ::augmentInterfacesVector({-4500,   50,  15600, 26500, 40500});
    auto velocities = ::augmentVelocityVector({       3500, 5900,  6400,  7500,  7900});
    auto slownesses = ::toSlownessVector(velocities);
    constexpr double sourceDepth{8000};
    constexpr double stationDepth{-1400};
    constexpr double stationOffset{20000};
    constexpr auto isAugmented{true};
    auto sourceLayer  = ::getLayer(sourceDepth,  interfaces, isAugmented);
    auto stationLayer = ::getLayer(stationDepth, interfaces, isAugmented);
    EXPECT_EQ(sourceLayer,  1); 
    EXPECT_EQ(stationLayer, 0); 

    double takeOffAngle{25};


    /// Snell's law
    double theta1Down = takeOffAngle*(M_PI/180);
    double theta2Down = std::asin((velocities[2]/velocities[1])*std::sin(theta1Down));
    double theta3Down = std::asin((velocities[3]/velocities[2])*std::sin(theta2Down));
    double theta3Up   = theta3Down;
    double theta2Up   = theta2Down;
    double theta1Up   = theta1Down; //M_PI - theta1Down;
    double theta0Up   = std::asin((velocities[0]/velocities[1]*std::sin(theta1Up)));
    //std::cout << theta0Up*180/M_PI << "," << theta1Up*180/M_PI << "," << theta1Down*180/M_PI << "," << theta0Up*180/M_PI << std::endl;
    //std::cout << theta0Up << std::endl;
    auto dz0Up   = interfaces[1] - stationDepth;
    auto dz1Down = interfaces[2] - sourceDepth;
    auto dz1Up   = interfaces[2] - interfaces[1];
    auto dz2Down = interfaces[3] - interfaces[2];
    auto dz2Up   = dz2Down;
    auto dz3Down = interfaces[4] - interfaces[3];
    auto dz3Up   = dz3Down;

    double dx1Down = dz1Down*std::tan(theta1Down); //10972.622618039548
    double dx2Down = dz2Down*std::tan(theta2Down); //8736.34515155134
    double dx3Down = dz3Down*std::tan(theta3Down); //15081.670759818773
    double dx3Up   = dz3Up*std::tan(theta3Up);
    double dx2Up   = dz2Up*std::tan(theta2Up);
    double dx1Up   = dz1Up*std::tan(theta1Up);
    double dx0Up   = dz0Up*std::tan(theta0Up);

    Path2D path1;
    {   
    Point2D point1{0, sourceDepth},
            point2{dx1Down, 15600.0},
            point3{dx1Down + dx1Up, 50},
            point4{dx1Down + dx1Up + dx0Up, stationDepth};
    Segment2D segment1, segment2, segment3, segment4;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[1]);
    segment1.setVelocityModelCellIndex(1);
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[1]);
    segment2.setVelocityModelCellIndex(1);
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[0]);
    segment3.setVelocityModelCellIndex(0);

    path1.set( {segment1, segment2, segment3} );
    }

    Path2D path2;
    {   
    Point2D point1{0, sourceDepth},
            point2{dx1Down, 15600.0},
            point3{dx1Down + dx2Down, 26500},
            point4{dx1Down + dx2Down + dx2Up, 15600},
            point5{dx1Down + dx2Down + dx2Up + dx1Up, 50},
            point6{dx1Down + dx2Down + dx2Up + dx1Up + dx0Up, stationDepth};
    Segment2D segment1, segment2, segment3, segment4, segment5;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[1]);
    segment1.setVelocityModelCellIndex(1);
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[2]);
    segment2.setVelocityModelCellIndex(2);
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[2]);
    segment3.setVelocityModelCellIndex(2);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[1]);
    segment4.setVelocityModelCellIndex(1);
    segment5.setStartAndEndPoint(std::pair {point5, point6});
    segment5.setSlowness(slownesses[0]);
    segment5.setVelocityModelCellIndex(0);

    path2.set( {segment1, segment2, segment3, segment4, segment5} );
    } 

    Path2D path3;
    {
    Point2D point1{0, sourceDepth},
            point2{dx1Down, 15600.0},
            point3{dx1Down + dx2Down, 26500},
            point4{dx1Down + dx2Down + dx3Down, 40500},
            point5{dx1Down + dx2Down + dx3Down + dx3Up, 26500},
            point6{dx1Down + dx2Down + dx3Down + dx3Up + dx2Up, 15600},
            point7{dx1Down + dx2Down + dx3Down + dx3Up + dx2Up + dx1Up, 50},
            point8{dx1Down + dx2Down + dx3Down + dx3Up + dx2Up + dx1Up + dx0Up, stationDepth};
    Segment2D segment1, segment2, segment3, segment4, segment5, segment6, segment7;
    segment1.setStartAndEndPoint(std::pair {point1, point2});
    segment1.setSlowness(slownesses[1]);
    segment1.setVelocityModelCellIndex(1);
    segment2.setStartAndEndPoint(std::pair {point2, point3});
    segment2.setSlowness(slownesses[2]);
    segment2.setVelocityModelCellIndex(2);
    segment3.setStartAndEndPoint(std::pair {point3, point4});
    segment3.setSlowness(slownesses[3]);
    segment3.setVelocityModelCellIndex(3);
    segment4.setStartAndEndPoint(std::pair {point4, point5});
    segment4.setSlowness(slownesses[3]);
    segment4.setVelocityModelCellIndex(3);
    segment5.setStartAndEndPoint(std::pair {point5, point6});
    segment5.setSlowness(slownesses[2]);
    segment5.setVelocityModelCellIndex(2);
    segment6.setStartAndEndPoint(std::pair {point6, point7});
    segment6.setSlowness(slownesses[1]);
    segment6.setVelocityModelCellIndex(1);
    segment7.setStartAndEndPoint(std::pair {point7, point8});
    segment7.setSlowness(slownesses[0]);
    segment7.setVelocityModelCellIndex(0);

    path3.set( {segment1, segment2, segment3, segment4, segment5, segment6, segment7} );
    }

    std::vector<Path2D> paths{path1, path2, path3};
    int iPath = 0;
    for (int layer = sourceLayer; layer < nLayers - 1; ++layer)
    {
        std::vector<::Segment> segments;
        auto returnCode = ::traceDownThenUp(interfaces,
                                            slownesses,
                                            takeOffAngle,
                                            sourceLayer,
                                            sourceDepth,
                                            stationLayer,
                                            stationOffset,
                                            stationDepth,
                                            layer,
                                            &segments,
                                            true,
                                            1);
        EXPECT_TRUE(returnCode == ReturnCode::Hit ||
                    returnCode == ReturnCode::UnderShot ||
                    returnCode == ReturnCode::OverShot);
        auto rayPath = ::toRayPath(segments);
        //for (const auto &segment : rayPath){std::cout << segment << std::endl;}
        //std::cout << std::endl;
        //for (const auto &segment : paths[iPath]){std::cout << segment << std::endl;}
        //std::cout << std::endl;
        //getchar();
        EXPECT_TRUE(rayPath == paths[iPath]);
        iPath = iPath + 1;
    }

}

TEST(Ray, LayerSolverHalfSpace)
{
    std::vector<double> interfaces{-4500};
    std::vector<double> velocities{3500};
    LayerSolver solver;
    solver.setVelocityModel(interfaces, velocities);
    solver.setSourceDepth(0);
    solver.setStationOffset(5000);
    solver.solve();
    EXPECT_TRUE(solver.haveRayPaths());
    auto rayPaths = solver.getRayPaths();
    EXPECT_EQ(rayPaths.size(), 1);

    Point2D point0{0, 0};
    Point2D point1{5000, -4500};
    Segment2D segment;
    segment.setStartAndEndPoint(std::pair {point0, point1}); 
    segment.setVelocity(velocities[0]);
    segment.setVelocityModelCellIndex(0);
    std::vector<Path2D> referenceRayPaths;
    Path2D path;
    path.open();
    path.append(segment);
    path.close();
    referenceRayPaths.push_back(std::move(path));
    EXPECT_TRUE(referenceRayPaths == rayPaths);
}


/*
TEST(Ray, LayerSolverVerticalUp)
{
    std::vector<double> interfaces{-4500, 40, 15600, 26500, 40500};
    std::vector<double> velocities{   3500, 5900,  6400,  7500,  7900};
    LayerSolver solver;
    solver.setVelocityModel(interfaces, velocities);
    solver.setSourceDepth(5000);
    solver.setStationOffset(0);
    solver.solve();
    auto rayPaths = solver.getRayPaths();
    EXPECT_EQ(rayPaths.size(), 1);

    // Ray goes up through first and second layer
    Path2D path;
    path.open();
    Point2D startPoint1{0, 5000};
    Point2D endPoint1{0, 40};
    Segment2D segment1;
    segment1.setStartAndEndPoint(std::pair {startPoint1, endPoint1});
    segment1.setVelocity(velocities[1]);
    segment1.setVelocityModelCellIndex(1);
    path.append(segment1);

    Point2D startPoint2{endPoint1};
    Point2D endPoint2{0, interfaces[0]};
    Segment2D segment2; 
    segment2.setStartAndEndPoint(std::pair {startPoint2, endPoint2});
    segment2.setVelocity(velocities[0]);
    segment2.setVelocityModelCellIndex(0);
    path.append(segment2);
    path.close();

    std::vector<Path2D> referenceRayPaths;
    referenceRayPaths.push_back(std::move(path));
    EXPECT_TRUE(referenceRayPaths == rayPaths);

    solver.setSourceDepth(40);
    solver.solve();
    auto dz1 = interfaces[1] - interfaces[0];
    EXPECT_EQ(solver.getRayPaths().at(0).size(), 1);
    auto travelTime = dz1/velocities[0];
    EXPECT_NEAR(travelTime, solver.getRayPaths().at(0).getTravelTime(), 1.e-10);

    // Deal with some edge edge cases
    solver.setSourceDepth(26500);
    solver.solve();
    auto dz3 = interfaces[3] - interfaces[2];
    auto dz2 = interfaces[2] - interfaces[1];
    travelTime = dz1/velocities[0] + dz2/velocities[1] + dz3/velocities[2];
    EXPECT_EQ(solver.getRayPaths().at(0).size(), 3);
    EXPECT_NEAR(travelTime, solver.getRayPaths().at(0).getTravelTime(), 1.e-10);
    //std::cout << solver.getRayPaths().at(0).size() << std::endl;
    //std::cout << travelTime << " " << solver.getRayPaths().at(0).getTravelTime() << std::endl;
    
}

TEST(Ray, LayerSolverVerticalDown)
{
    std::vector<double> interfaces{-4500, 40, 15600, 26500, 40500};
    std::vector<double> velocities{   3500, 5900,  6400,  7500,  7900};
    LayerSolver solver;
    solver.setVelocityModel(interfaces, velocities);
    solver.setSourceDepth(-4000);
    solver.setStationOffsetAndDepth(0, 16000);
    solver.solve();
    auto rayPaths = solver.getRayPaths();
    EXPECT_EQ(rayPaths.size(), 1);

    // Ray goes up through first and second layer
    Path2D path;
    path.open();
    Point2D startPoint1{0, -4000};
    Point2D endPoint1{0, interfaces[1]};
    Segment2D segment1;
    segment1.setStartAndEndPoint(std::pair {startPoint1, endPoint1});
    segment1.setVelocity(velocities[0]);
    segment1.setVelocityModelCellIndex(0);
    path.append(segment1);

    Point2D startPoint2{endPoint1};
    Point2D endPoint2{0, interfaces[2]};
    Segment2D segment2;
    segment2.setStartAndEndPoint(std::pair {startPoint2, endPoint2});
    segment2.setVelocity(velocities[1]);
    segment2.setVelocityModelCellIndex(1);
    path.append(segment2);

    Point2D startPoint3{endPoint2};
    Point2D endPoint3{0, 16000};
    Segment2D segment3;
    segment3.setStartAndEndPoint(std::pair {startPoint3, endPoint3});
    segment3.setVelocity(velocities[2]);
    segment3.setVelocityModelCellIndex(2);
    path.append(segment3);
    path.close();

    std::vector<Path2D> referenceRayPaths;
    referenceRayPaths.push_back(std::move(path));
    EXPECT_TRUE(referenceRayPaths == rayPaths);
}
*/

TEST(Ray, FirstArrivalSourceStationSameDepth)
{
    constexpr double depth{-2000};
    std::vector<double> stationOffsets{1000,
                                       5000,
                                       7500,
                                       10000,
                                       15000,
                                       20000,
                                       30000,
                                       40000,
                                       50000,
                                       60000,
                                       80000,
                                       100000,
                                       120000};
    // TODO travel times from numerical ray paths seem more accurate.
    // Why is eikonal fast by 0.01 seconds?
    std::vector<double> referenceTravelTimes{0.285714285714,
                                             1.42857142857,
                                             2.14285714286,
                                             2.62783136508,
                                             3.4752889922,
                                             4.32274661932,
                                             6.01766187355,
                                             7.71257712779,
                                             9.40749238203,
                                             11.1024076363,
                                             14.4922381447,
                                             17.8820686532,
                                             21.2718991617};
    std::vector<double> interfaces{-4500, 50, 15600, 26500, 40500};
    std::vector<double> velocities{   3500, 5900,  6400,  7500,  7900};
//    ::createReferenceSolution(depth, depth,
//                              interfaces, velocities, stationOffsets);
    LayerSolver solver;
    solver.setVelocityModel(interfaces, velocities);
    solver.setSourceDepth(depth);
    for (int iOffset = 0;
         iOffset < static_cast<int> (stationOffsets.size());
         ++iOffset)
    {
        solver.setStationOffsetAndDepth(stationOffsets[iOffset], depth);
        solver.solve();
        auto rayPaths = solver.getRayPaths();
        EXPECT_NEAR(referenceTravelTimes.at(iOffset),
                    rayPaths.at(0).getTravelTime(),
                    0.02); 
        std::cout << "Eikonal vs computed travel time: "
                  << referenceTravelTimes.at(iOffset) << "," 
                  << rayPaths.at(0).getTravelTime() << ","
                  << rayPaths.at(0).getTakeOffAngle() << std::endl;
        //for (const auto &segment : rayPaths[0])
        //{
        //    std::cout << segment << std::endl;
        //}
    }
}

TEST(Ray, FirstArrivalSourceDeeperThanStation)
{
    constexpr double stationDepth{-1400};
    constexpr double sourceDepth{18000};
    std::vector<double> stationOffsets{0,
                                       1000,
                                       5000,
                                       7500,
                                       10000,
                                       15000,
                                       20000,
                                       30000,
                                       40000,
                                       50000,
                                       60000,
                                       80000,
                                       100000,
                                       120000};
    // TODO travel times from numerical ray paths seem more accurate.
    // Why is eikonal fast by 0.01 seconds?
    std::vector<double> referenceTravelTimes{3.42088702324,
                                             3.42541195061,
                                             3.53043481827,
                                             3.66244781813,
                                             3.83923874418,
                                             4.30245895143,
                                             4.87382792766,
                                             6.20981834293,
                                             7.67981680992,
                                             9.20774986021,
                                             10.7572145407,
                                             13.8730672134,
                                             16.8996388956,
                                             19.5663055623};
    std::vector<double> interfaces{-4500, 50, 15600, 26500, 40500};
    std::vector<double> velocities{   3500, 5900,  6400,  7500,  7900};
//    ::createReferenceSolution(sourceDepth, stationDepth,
//                              interfaces, velocities, stationOffsets);
    LayerSolver solver;
    solver.setVelocityModel(interfaces, velocities);
    solver.setSourceDepth(sourceDepth);
    for (int iOffset = 0;
         iOffset < static_cast<int> (stationOffsets.size());
         ++iOffset)
    {
        solver.setStationOffsetAndDepth(stationOffsets[iOffset], stationDepth);
        solver.solve();
        auto rayPaths = solver.getRayPaths();
        EXPECT_NEAR(referenceTravelTimes.at(iOffset),
                    rayPaths.at(0).getTravelTime(),
                    0.02); 
        std::cout << "Eikonal vs computed travel time for test 2: "
                  << referenceTravelTimes.at(iOffset) << ","
                  << rayPaths.at(0).getTravelTime() << ","
                  << rayPaths.at(0).getTakeOffAngle() << std::endl;
    }
}

}

