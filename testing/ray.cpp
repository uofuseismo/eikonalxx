#include <iostream>
#include <cmath>
#include <vector>
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/point3d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/segment3d.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Ray, Point2D)
{
    Ray::Point2D point;
    const double x{5};
    const double z{7};

    point.setPositionInX(x);
    point.setPositionInZ(z);

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
    Ray::Point3D point;
    const double x{5};
    const double y{6};
    const double z{7};

    point.setPositionInX(x);
    point.setPositionInY(y);
    point.setPositionInZ(z);

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
    EXPECT_NEAR(copy.getVelocity(), velocity, 1.e-15);
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
    EXPECT_NEAR(copy.getVelocity(), velocity, 1.e-15);
    EXPECT_NEAR(copy.getSlowness(), 1./velocity, 1.e-15);
    EXPECT_NEAR(copy.getLength(),   length, 1.e-15);
    EXPECT_NEAR(copy.getTravelTime(), travelTime, 1.e-15);
    EXPECT_EQ(copy.getVelocityModelCellIndex(), cellIndex);

    segment.clear();
    EXPECT_FALSE(segment.haveStartAndEndPoint());
    EXPECT_FALSE(segment.haveVelocity());
    EXPECT_FALSE(segment.haveVelocityModelCellIndex());
}

}
