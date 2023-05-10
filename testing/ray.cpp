#include <iostream>
#include <fstream>
#include <vector>
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/point3d.hpp"
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

}
