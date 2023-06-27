#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/layerSolver.hpp"
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

}

