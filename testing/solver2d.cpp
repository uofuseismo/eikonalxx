#include <iostream>
#include <fstream>
#include <vector>
#include "eikonalxx/solver2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "private/solverUtilities2d.hpp"
#include "private/solver2d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Solver2D, numberOfLevels)
{
    int nx = 5;
    int nz = 7;
    int nLevels = nx + nz - 1;
    EXPECT_EQ(computeNumberOfLevels(nx, nz), nLevels);
}

TEST(Solver2D, sweepSigns)
{
    int ixShift, izShift, signX, signZ;
    getSweepFiniteDifferenceSigns<SweepNumber2D::SWEEP1>(&ixShift, &izShift,
                                                         &signX, &signZ);
    EXPECT_EQ(ixShift,-1);
    EXPECT_EQ(izShift,-1);
    EXPECT_EQ(signX, 1);
    EXPECT_EQ(signZ, 1);

    getSweepFiniteDifferenceSigns<SweepNumber2D::SWEEP2>(&ixShift, &izShift,
                                                         &signX, &signZ);
    EXPECT_EQ(ixShift, 1);
    EXPECT_EQ(izShift,-1);
    EXPECT_EQ(signX,-1);
    EXPECT_EQ(signZ, 1);

    getSweepFiniteDifferenceSigns<SweepNumber2D::SWEEP3>(&ixShift, &izShift,
                                                         &signX, &signZ);
    EXPECT_EQ(ixShift,-1);
    EXPECT_EQ(izShift, 1);
    EXPECT_EQ(signX, 1);
    EXPECT_EQ(signZ,-1);

    getSweepFiniteDifferenceSigns<SweepNumber2D::SWEEP4>(&ixShift, &izShift,
                                                         &signX, &signZ);
    EXPECT_EQ(ixShift, 1);
    EXPECT_EQ(izShift, 1);
    EXPECT_EQ(signX,-1);
    EXPECT_EQ(signZ,-1);
}

//----------------------------------------------------------------------------//

TEST(Solver2D, analyticSolutions)
{
    double dx = 100.;
    double dz = 200.;
    double xSrc = 50;
    double zSrc = 100;
    double slowness = 1./5000.;
    double t, dtdx, dtdz;
    
    double distance = std::sqrt(xSrc*xSrc + zSrc*zSrc);
    double tref = distance*slowness;
    auto theta = std::atan2(dz/2, dx/2);
    // dT/dx and dT/dz are really just apparent slownesses
    auto sx = std::cos(theta)*slowness;
    auto sz = std::sin(theta)*slowness; 

    // Up and left
    int ix = 0;
    int iz = 0;
    t = computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness);
    EXPECT_NEAR(t, tref, 1.e-10);
//std::cout << tref << std::endl;
    computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness,
                                &t, &dtdx, &dtdz); 
    EXPECT_NEAR(t, tref, 1.e-10);
    EXPECT_NEAR(-sx, dtdx, 1.e-10);
    EXPECT_NEAR(-sz, dtdz, 1.e-10);
//std::cout << t << std::endl;
//std::cout << dtdx << std::endl;
//std::cout << dtdz << std::endl;

    // Up and right
    ix = 1;
    iz = 0;
    t = computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness);
    EXPECT_NEAR(t, tref, 1.e-10);
    computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness,
                                &t, &dtdx, &dtdz);
    EXPECT_NEAR(t, tref, 1.e-10);
    EXPECT_NEAR( sx, dtdx, 1.e-10);
    EXPECT_NEAR(-sz, dtdz, 1.e-10);

    // Down and right
    ix = 1;
    iz = 1;
    t = computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness);
    EXPECT_NEAR(t, tref, 1.e-10);
    computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness,
                                &t, &dtdx, &dtdz);
    EXPECT_NEAR(t, tref, 1.e-10);
    EXPECT_NEAR( sx, dtdx, 1.e-10);
    EXPECT_NEAR( sz, dtdz, 1.e-10);

    // Down and left
    ix = 0;
    iz = 1;
    t = computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness);
    EXPECT_NEAR(t, tref, 1.e-10);
    computeAnalyticalTravelTime(ix, iz, dx, dz, xSrc, zSrc, slowness,
                                &t, &dtdx, &dtdz);
    EXPECT_NEAR(t, tref, 1.e-10);
    EXPECT_NEAR(-sx, dtdx, 1.e-10);
    EXPECT_NEAR( sz, dtdz, 1.e-10);

    // Edge case
    computeAnalyticalTravelTime(0, 0, 1., 1., 0., 0., slowness,
                                &t, &dtdx, &dtdz);
    EXPECT_NEAR(t, 0, 1.e-10);
    EXPECT_NEAR(dtdx, 0, 1.e-14);
    EXPECT_NEAR(dtdz, 0, 1.e-14);
    //EXPECT_TRUE(std::isnan(dtdx));
    //EXPECT_TRUE(std::isnan(dtdz));
}

//----------------------------------------------------------------------------//

/*
TEST(Solver2D, testFiniteDifferences)
{
    const double huge = 1.e10;
    const double defaultVel = 4000;
    double dx, dz = 100;
    double tupd = 0;
    double t0, t1, t2, t3 = 0;
    double s0, s1, s2, s3 = 1./defaultVel;
    double dx_dz = dx/dz;
    double dz_dx = dz/dx;
    double cosTheta = dx/std::sqrt(dx*dx + dz*dz);
    double sinTheta = dz/std::sqrt(dx*dx + dz*dz);
    int sphericalRadius = 5;
    // Test critical refraction in x and z
    dx = 100;
    dz = 200;
    t1 = 1;
    t2 = huge; // deactivate 3 point operators
    t3 = huge;
    s0 = s1 = s2 = s3 = 1./defaultVel;
    s3 = 1/5000.;
    cosTheta = dx/std::sqrt(dx*dx + dz*dz);
    sinTheta = dz/std::sqrt(dx*dx + dz*dz);
    dx_dz = dx/dz;
    dz_dx = dz/dx;
    //tupd = cartesianFiniteDifferenceTwoPoint(dx, dz, s0, s1, s3, t1, t3);
    tupd = cartesianFiniteDifference(huge, 
                                     dx, dz,
                                     dx_dz, dz_dx,
                                     cosTheta, sinTheta,
                                     s0, s1, s3,
                                     t1, t2, t3); 
std::cout << "Traveltime update in x: " << tupd << std::endl;
    EXPECT_NEAR(tupd, t1 + dx*s3, 1.e-10); 

    t1 = huge;
    t2 = huge; // deactivate 3 point operators
    t3 = 2;
    s0 = s1 = s2 = s3 = 1./defaultVel;
    s1 = 1/4500.;
    //tupd = cartesianFiniteDifferenceTwoPoint(dx, dz, s0, s1, s3, t1, t3); 
    tupd = cartesianFiniteDifference(huge,
                                     dx, dz,
                                     dx_dz, dz_dx,
                                     cosTheta, sinTheta,
                                     s0, s1, s3,
                                     t1, t2, t3);
std::cout << "Traveltime update in z: " << tupd << std::endl;
    EXPECT_NEAR(tupd, t3 + dz*s1, 1.e-10); 

    s0 = s1 = s2 = s3 = 1./defaultVel;

    // Podvin and Lecomte
    dx = 100;
    dz = 109;
    s0 = 1./defaultVel;
    t2 = 1;
    t1 = t2 + s0*std::sqrt(dx*dx + dz*dz)/2.;
    t3 = t1;  
    s1 = s2 = s3 = 1/3000.; // Deactivate critical refraction
    cosTheta = dx/std::sqrt(dx*dx + dz*dz);
    sinTheta = dz/std::sqrt(dx*dx + dz*dz); 
    dx_dz = dx/dz;
    dz_dx = dz/dx;
    tupd = cartesianFiniteDifference(huge,
                                     dx, dz,
                                     dx_dz, dz_dx,
                                     cosTheta, sinTheta,
                                     s0, s1, s3,
                                     t1, t2, t3); 
std::cout << tupd - 1.03684367300601 << tupd << " " << 1 + s0*std::sqrt(dx*dx + dz*dz) << std::endl;
//tupd = cartesianFiniteDifference(huge, dx, s0, s1, s3, t1, t2, t3);   
//std::cout << tupd << std::endl;

    dx = 100;
    dz = 100;
    SweepNumber sweep = SweepNumber2D::Sweep1;
    double sourceSlowness = 1./defaultVel;
    double xSourceOffset = 0;
    double zSourceOffset = 0;
    auto iSrcX = static_cast<int> (xSourceOffset/dx);
    auto iSrcZ = static_cast<int> (zSourceOffset/dz);
    int ix = 1;
    int iz = 1;
    s0 = s1 = s2 = s3 = 1./defaultVel; 
    t2 = 0;
    t1 = ix*dz*s0;
    t3 = iz*dx*s0; 
    int ixShift, izShift, signX, signZ;
    getSweepFiniteDifferenceSigns(SweepNumber2D::Sweep1,
                                  &ixShift, &izShift,
                                  &signX, &signZ);
    tupd = sphericalFiniteDifference(sphericalRadius,
                                     huge, dx, sourceSlowness,
                                     ix, iz,
                                     ixShift, izShift,
                                     signX, signZ,
                                     iSrcX, iSrcZ,
                                     xSourceOffset, zSourceOffset,
                                     s0, s1, s3,
                                     t1, t2, t3);
    std::cout << tupd << " " << t2 + s0*std::sqrt(dx*dx + dz*dz) << std::endl;
}
*/

//----------------------------------------------------------------------------//

TEST(Solver2D, levelSetIndices)
{
    int nx = 3;
    int nz = 5;
    auto nLevels = computeNumberOfLevels(nx, nz);
    auto levelOffset = makeLevelOffset(nx, nz);
    EXPECT_EQ(levelOffset[0], 0);
    EXPECT_EQ(static_cast<int> (levelOffset.size()), nLevels + 1);
    std::vector<std::pair<int, int>> refLevels
       { {0, 0}, // (ix,iz)
         {1, 0}, {0, 1}, 
         {2, 0}, {1, 1}, {0, 2},
         {2, 1}, {1, 2}, {0, 3},
         {2, 2}, {1, 3}, {0, 4},
         {2, 3}, {1, 4},
         {2, 4} };
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        int j = 0;
        for (int level=0; level<nLevels; ++level)
        {
            int i0, i1;
            int ix, iz;
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx=i0; indx<i1; ++indx)
            {
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
                else
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
            //std::cout << "(level,ix,iz)" << level << " " << ix << " " << iz << " " << refLevels[j].first <<  "," <<  refLevels[j].second << std::endl;
                j = j + 1;
            }
            EXPECT_EQ(levelOffset[level+1], j);
        }
        EXPECT_EQ(j, static_cast<int> (refLevels.size()));
    }
    // Transpose previous case
    nx = 5;
    nz = 3;
    nLevels = computeNumberOfLevels(nx, nz);
    levelOffset = makeLevelOffset(nx, nz);
    EXPECT_EQ(levelOffset[0], 0);
    EXPECT_EQ(static_cast<int> (levelOffset.size()), nLevels + 1);
    refLevels = { {0, 0}, // (ix,iz)
                  {1, 0}, {0, 1},
                  {2, 0}, {1, 1}, {0, 2},
                  {3, 0}, {2, 1}, {1, 2},
                  {4, 0}, {3, 1}, {2, 2},
                  {4, 1}, {3, 2},
                  {4, 2} };
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        int j = 0;
        for (int level=0; level<nLevels; ++level)
        {
            int i0, i1;
            int ix, iz;
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx=i0; indx<i1; ++indx)
            {
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
                else
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
                j = j + 1;
            }
            EXPECT_EQ(levelOffset[level+1], j);
        }
        EXPECT_EQ(j, static_cast<int> (refLevels.size()));
    }

    // Do even grid-size case
    nx = 4;
    nz = 4;
    nLevels = computeNumberOfLevels(nx, nz);
    levelOffset = makeLevelOffset(nx, nz);
    EXPECT_EQ(levelOffset[0], 0);
    EXPECT_EQ(static_cast<int> (levelOffset.size()), nLevels + 1);
    refLevels = { {0, 0},
                  {1, 0}, {0, 1},
                  {2, 0}, {1, 1}, {0, 2},
                  {3, 0}, {2, 1}, {1, 2}, {0, 3},
                  {3, 1}, {2, 2}, {1, 3},
                  {3, 2}, {2, 3}, 
                  {3, 3} };
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        int j = 0;
        for (int level = 0; level < nLevels; ++level)
        {
            int i0, i1;
            int ix, iz;
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx = i0; indx < i1; ++indx)
            {
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, refLevels[j].second);
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
                else
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    EXPECT_EQ(ix, nx - 1 - refLevels[j].first);
                    EXPECT_EQ(iz, nz - 1 - refLevels[j].second);
                }
                j = j + 1;
            }
            EXPECT_EQ(levelOffset[level+1], j);
        }
        EXPECT_EQ(j, static_cast<int> (refLevels.size()));
    }
}

//============================================================================//

TEST(Solver2D, gridSweepToLevelIndex)
{
    // Equal grid size
    int nx = 12;
    int nz = 12;
    auto nLevels = computeNumberOfLevels(nx, nz);
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        for (int level=0; level<nLevels; ++level)
        {
            int i0, i1;
            int ix, iz;
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx = i0; indx < i1; ++indx)
            {
                int levelTest, indxTest =-1;
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP1>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP2>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP3>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP4)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP4>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                EXPECT_NE(levelTest, -1);
                EXPECT_NE(indxTest, -1);
                EXPECT_EQ(level, levelTest);
                EXPECT_EQ(indx, indxTest); 
             }
        }
    }
    // Unequal grid size
    nx = 52;
    nz = 71;
    nLevels = computeNumberOfLevels(nx, nz);
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        for (int level = 0; level < nLevels; ++level)
        {
            int i0, i1;
            int ix, iz;
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx = i0; indx < i1; ++indx)
            {
                int levelTest, indxTest =-1;
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP1>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP2>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP3>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                else if (sweep == SweepNumber2D::SWEEP4)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(level, indx,
                                                                 nx, nz,
                                                                 &ix, &iz);
                    gridSweepToLevelIndex<SweepNumber2D::SWEEP4>(ix, iz, nx, nz,
                                                                 &levelTest,
                                                                 &indxTest);
                }
                EXPECT_NE(levelTest, -1);
                EXPECT_NE(indxTest, -1);
                EXPECT_EQ(level, levelTest);
                EXPECT_EQ(indx, indxTest);
             }
        }
    }
}

//----------------------------------------------------------------------------//

TEST(Solver2D, sweepToGridTravelTimeIndices)
{
    int nx = 5;
    int nz = 3;
    auto nLevels = computeNumberOfLevels(nx, nz);
    auto levelOffset = makeLevelOffset(nx, nz);    
    std::vector<double> tTimes{1,  2, 3, 4, 5,
                               6,  7, 8, 9, 10,
                               11,12,13,14,15};
    std::vector<double> tRef0{ 1, 1, 1, 1,
                               2, 1, 1, 2,   6, 6, 1, 1,
                               3, 2, 2, 3,   7, 6, 1, 2,  11,11, 6, 6,
                               4, 3, 3, 4,   8, 7, 2, 3,  12,11, 6, 7,
                               5, 4, 4, 5,   9, 8, 3, 4,  13,12, 7, 8,
                              10, 9, 4, 5,  14,13, 8, 9,
                              15,14, 9,10};
    std::vector<double> tRef1{ 5, 5, 5, 5,
                               4, 5, 5, 4,  10,10, 5, 5,
                               3, 4, 4, 3,   9,10, 5, 4,  15,15,10,10,
                               2, 3, 3, 2,   8, 9, 4, 3,  14,15,10, 9,
                               1, 2, 2, 1,   7, 8, 3, 2,  13,14, 9, 8,
                               6, 7, 2, 1,  12,13, 8, 7,
                              11,12, 7, 6};
    std::vector<double> tRef2{11,11,11,11,
                              12,11,11,12,   6, 6,11,11,
                              13,12,12,13,   7, 6,11,12,  1, 1, 6, 6,
                              14,13,13,14,   8, 7,12,13,  2, 1, 6, 7,
                              15,14,14,15,   9, 8,13,14,  3, 2, 7, 8,
                              10, 9,14,15,   4, 3, 8, 9,
                               5, 4, 9,10}; 
    std::vector<double> tRef3{15,15,15,15,
                              14,15,15,14,  10,10,15,15,
                              13,14,14,13,   9,10,15,14,  5, 5,10,10,
                              12,13,13,12,   8, 9,14,13,  4, 5,10, 9,
                              11,12,12,11,   7, 8,13,12,  3, 4, 9, 8,
                               6, 7,12,11,   2, 3, 8, 7,
                               1, 2, 7, 6};
    int j = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        int iStart, iEnd;
        int it0, it1, it2, it3;
        getLevelStartStopIndices(nx, nz, level, &iStart, &iEnd);
        for (int indx = iStart; indx < iEnd; ++indx)
        {
            sweepLevelIndexToTravelTimeIndices<SweepNumber2D::SWEEP1>(
                                               level, indx,
                                               nx, nz, 
                                               &it0, &it1, &it2, &it3);
            auto jt0 = sweepLevelIndexToIndex<SweepNumber2D::SWEEP1>(level,
                                                                     indx,
                                                                     nx, nz);
            EXPECT_NEAR(tTimes[it0], tRef0.at(j),   1.e-10);
            EXPECT_NEAR(tTimes[it1], tRef0.at(j+1), 1.e-10);
            EXPECT_NEAR(tTimes[it2], tRef0.at(j+2), 1.e-10);
            EXPECT_NEAR(tTimes[it3], tRef0.at(j+3), 1.e-10);
            EXPECT_NEAR(tTimes[jt0], tRef0.at(j),   1.e-10);
//     std::cout << tTimes.at(it0) << " " << tRef0.at(j) << std::endl; /// << " " << it1 << " " << it2 << " " << it3 << "|"; 
            j = j + 4;
        }
//std::cout << std::endl;
    }
    EXPECT_EQ(j, static_cast<int> (tRef0.size()));

    j = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        int iStart, iEnd;
        int it0, it1, it2, it3;
        getLevelStartStopIndices(nx, nz, level, &iStart, &iEnd);
        for (int indx = iStart; indx < iEnd; ++indx)
        {
            sweepLevelIndexToTravelTimeIndices<SweepNumber2D::SWEEP2>(
                                               level, indx,
                                               nx, nz,
                                               &it0, &it1, &it2, &it3);
            auto jt0 = sweepLevelIndexToIndex<SweepNumber2D::SWEEP2>(level,
                                                                     indx,
                                                                     nx, nz);
            EXPECT_NEAR(tTimes[it0], tRef1.at(j),   1.e-10);
            EXPECT_NEAR(tTimes[it1], tRef1.at(j+1), 1.e-10);
            EXPECT_NEAR(tTimes[it2], tRef1.at(j+2), 1.e-10);
            EXPECT_NEAR(tTimes[it3], tRef1.at(j+3), 1.e-10);
            EXPECT_NEAR(tTimes[jt0], tRef1.at(j),   1.e-10);
            j = j + 4;
        }
    }
    EXPECT_EQ(j, static_cast<int> (tRef1.size()));

    j = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        int iStart, iEnd;
        int it0, it1, it2, it3;
        getLevelStartStopIndices(nx, nz, level, &iStart, &iEnd);
//std::cout << level << std::endl;
        for (int indx = iStart; indx < iEnd; ++indx)
        {
            sweepLevelIndexToTravelTimeIndices<SweepNumber2D::SWEEP3>(
                                               level, indx,
                                               nx, nz,
                                               &it0, &it1, &it2, &it3);
            auto jt0 = sweepLevelIndexToIndex<SweepNumber2D::SWEEP3>(level,
                                                                     indx,
                                                                     nx, nz);
//std::cout << it0 << " " << it1 << " " << it2 << " " << it3 << "|";
            EXPECT_NEAR(tTimes[it0], tRef2.at(j),   1.e-10);
            EXPECT_NEAR(tTimes[it1], tRef2.at(j+1), 1.e-10);
            EXPECT_NEAR(tTimes[it2], tRef2.at(j+2), 1.e-10);
            EXPECT_NEAR(tTimes[it3], tRef2.at(j+3), 1.e-10);
            EXPECT_NEAR(tTimes[jt0], tRef2.at(j),   1.e-10);
            j = j + 4;
        }
//std::cout << std::endl;
    }
    EXPECT_EQ(j, static_cast<int> (tRef2.size()));

    j = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        int iStart, iEnd;
        int it0, it1, it2, it3;
        getLevelStartStopIndices(nx, nz, level, &iStart, &iEnd);
//std::cout << level << std::endl;
        for (int indx = iStart; indx < iEnd; ++indx)
        {
            sweepLevelIndexToTravelTimeIndices<SweepNumber2D::SWEEP4>(
                                               level, indx,
                                               nx, nz,
                                               &it0, &it1, &it2, &it3);
            auto jt0 = sweepLevelIndexToIndex<SweepNumber2D::SWEEP4>(level,
                                                                     indx,
                                                                     nx, nz);
//std::cout << it0 << " " << it1 << " " << it2 << " " << it3 << "|";
            EXPECT_NEAR(tTimes[it0], tRef3.at(j),   1.e-10);
            EXPECT_NEAR(tTimes[it1], tRef3.at(j+1), 1.e-10);
            EXPECT_NEAR(tTimes[it2], tRef3.at(j+2), 1.e-10);
            EXPECT_NEAR(tTimes[it3], tRef3.at(j+3), 1.e-10);
            EXPECT_NEAR(tTimes[jt0], tRef3.at(j),   1.e-10);
            j = j + 4;
        }
//std::cout << std::endl;
    }
    EXPECT_EQ(j, static_cast<int> (tRef3.size()));
} 

//----------------------------------------------------------------------------//

TEST(Solver2D, setVelocityModel)
{
    int nx = 5;
    int nz = 4;
    int ncx = nx - 1;
    int ncz = nz - 1;
    std::vector<double> slow{1, 2,  3,  4, 
                             5, 6,  7,  8,
                             9, 10, 11, 12};
    std::vector<double> slowRef{1, 1, 1, 1};
    //std::vector<double> sweepSlowness(4*ngrid, 0);
    std::vector<SweepSlowness2D<double>> sweepSlowness;
    std::vector<double> sweepSlowRef1{ 1, 1, 1, 1,
                                       1, 2, 2, 1,   1, 1, 5, 5,
                                       2, 3, 3, 2,   1, 2, 6, 5,   5, 5, 9, 9,
                                       3, 4, 4, 3,   2, 3, 7, 6,   5, 6,10, 9,  9, 9, 9, 9,
                                       4, 4, 4, 4,   3, 4, 8, 7,   6, 7,11,10,  9,10,10, 9,
                                       4, 4, 8, 8,   7, 8,12,11,  10,11,11,10,
                                       8, 8,12,12,  11,12,12,11, 
                                      12,12,12,12};
                                      

    std::vector<double> sweepSlowRef2{ 4, 4, 4, 4,
                                       4, 3, 3, 4,  4, 4, 8, 8,
                                       3, 2, 2, 3,  4, 3, 7, 8,  8, 8,12,12,
                                       2, 1, 1, 2,  3, 2, 6, 7,  8, 7,11,12,  12,12,12,12,
                                       1, 1, 1, 1,  2, 1, 5, 6,  7, 6,10,11,  12,11,11,12,
                                       1, 1, 5, 5,  6, 5, 9,10, 11,10,10,11,
                                       5, 5, 9, 9, 10, 9, 9,10,
                                       9, 9, 9, 9};
    
    std::vector<double> sweepSlowRef3{ 9, 9, 9, 9,
                                       9,10,10, 9,  9, 9, 5, 5,
                                      10,11,11,10,  9,10, 6, 5,  5, 5, 1, 1,
                                      11,12,12,11, 10,11, 7, 6,  5, 6, 2, 1,  1, 1, 1, 1,
                                      12,12,12,12, 11,12, 8, 7,  6, 7, 3, 2,  1, 2, 2, 1,
                                      12,12, 8, 8,  7, 8, 4, 3,  2, 3, 3, 2, 
                                       8, 8, 4, 4,  3, 4, 4, 3,
                                       4, 4, 4, 4};

    std::vector<double> sweepSlowRef4{12,12,12,12,
                                      12,11,11,12, 12,12, 8, 8,
                                      11,10,10,11, 12,11, 7, 8,  8, 8, 4, 4,
                                      10, 9, 9,10, 11,10, 6, 7,  8, 7, 3, 4,  4, 4, 4, 4,
                                       9, 9, 9, 9, 10, 9, 5, 6,  7, 6, 2, 3,  4, 3, 3, 4,
                                       9, 9, 5, 5,  6, 5, 1, 2,  3, 2, 2, 3,
                                       5, 5, 1, 1,  2, 1, 1, 2,
                                       1, 1, 1, 1};
                                      
    //auto sweep = SweepNumber2D::SWEEP1;
    auto nLevels = computeNumberOfLevels(nx, nz);
    auto levelOffset = makeLevelOffset(nx, nz);
    sweepSlowness.resize(nLevels);
    for (int level=0; level<nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        sweepSlowness[level].allocate(nNodes);
    }
    slownessToSweepSlowness<double, SweepNumber2D::SWEEP1>(nLevels,
                                                   nx, nz, ncx, ncz,
                                                   slow.data(),
                                                   sweepSlowness.data());
    int i0 = 0;
    for (int level=0; level<nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        for (int i = 0; i < nNodes; ++i)
        {
            EXPECT_NEAR(sweepSlowness[level].s0[i],
                        sweepSlowRef1.at(i0),   1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s1[i],
                        sweepSlowRef1.at(i0+1), 1.e-10);
            //EXPECT_NEAR(sweepSlowness[level].s2[i],
            //            sweepSlowRef1.at(i0+2), 1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s3[i],
                        sweepSlowRef1.at(i0+3), 1.e-10);
            i0 = i0 + 4;
        }
    }
    EXPECT_EQ(i0, static_cast<int> (sweepSlowRef1.size()));
    
    //sweep = SweepNumber2D::SWEEP2;
    slownessToSweepSlowness<double, SweepNumber2D::SWEEP2>(nLevels,
                                                   nx, nz, ncx, ncz, 
                                                   slow.data(),
                                                   sweepSlowness.data());
    i0 = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        for (int i = 0; i < nNodes; ++i)
        {
            EXPECT_NEAR(sweepSlowness[level].s0[i],
                        sweepSlowRef2.at(i0),   1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s1[i],
                        sweepSlowRef2.at(i0+1), 1.e-10);
            //EXPECT_NEAR(sweepSlowness[level].s2[i],
            //            sweepSlowRef2.at(i0+2), 1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s3[i],
                        sweepSlowRef2.at(i0+3), 1.e-10);
            i0 = i0 + 4;
        }
    }
    EXPECT_EQ(i0, static_cast<int> (sweepSlowRef2.size()));

    //sweep = SweepNumber2D::SWEEP3;
    slownessToSweepSlowness<double, SweepNumber2D::SWEEP3>(nLevels,
                                                   nx, nz, ncx, ncz, 
                                                   slow.data(),
                                                   sweepSlowness.data());
    i0 = 0;
    for (int level = 0; level < nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        for (int i = 0; i < nNodes; ++i)
        {
            EXPECT_NEAR(sweepSlowness[level].s0[i],
                        sweepSlowRef3.at(i0),   1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s1[i],
                        sweepSlowRef3.at(i0+1), 1.e-10);
            //EXPECT_NEAR(sweepSlowness[level].s2[i],
            //            sweepSlowRef3.at(i0+2), 1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s3[i],
                        sweepSlowRef3.at(i0+3), 1.e-10);
            i0 = i0 + 4;
        }
    }
    EXPECT_EQ(i0, static_cast<int> (sweepSlowRef3.size()));

    //sweep = SweepNumber2D::SWEEP4;
    slownessToSweepSlowness<double, SweepNumber2D::SWEEP4>(nLevels,
                                                   nx, nz, ncx, ncz, 
                                                   slow.data(),
                                                   sweepSlowness.data());
    i0 = 0;
    for (int level=0; level<nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        for (int i = 0; i < nNodes; ++i)
        {
            EXPECT_NEAR(sweepSlowness[level].s0[i],
                        sweepSlowRef4.at(i0),   1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s1[i],
                        sweepSlowRef4.at(i0+1), 1.e-10);
            //EXPECT_NEAR(sweepSlowness[level].s2[i],
            //            sweepSlowRef4.at(i0+2), 1.e-10);
            EXPECT_NEAR(sweepSlowness[level].s3[i],
                        sweepSlowRef4.at(i0+3), 1.e-10);
            i0 = i0 + 4;
        }
    }
    EXPECT_EQ(i0, static_cast<int> (sweepSlowRef4.size()));

    // Do a more involved test 
    //sweep = SweepNumber2D::SWEEP1;
    nx = 242;
    nz = 323;
    ncx = nx - 1;
    ncz = nz - 1;
    sweepSlowness.clear(); //resize(4*ngrid, 0);
    nLevels = computeNumberOfLevels(nx, nz);
    levelOffset = makeLevelOffset(nx, nz);
    sweepSlowness.resize(nLevels);
    for (int level = 0; level < nLevels; ++level)
    {
        auto nNodes = levelOffset[level+1] - levelOffset[level];
        sweepSlowness[level].allocate(nNodes);
    }
    slow.resize(ncx*ncz, 0);
    for (int iz = 0; iz < ncz; ++iz)
    {
        for (int ix=0; ix<ncx; ++ix)
        {
            slow[iz*ncx + ix] = iz*ncx + ix + 1;
        }
    }
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        if (is == 0)
        {
            slownessToSweepSlowness<double, SweepNumber2D::SWEEP1>(nLevels,
                                                          nx, nz, ncx, ncz,
                                                          slow.data(),
                                                          sweepSlowness.data());
        }
        else if (is == 1)
        {
            slownessToSweepSlowness<double, SweepNumber2D::SWEEP2>(nLevels, 
                                                          nx, nz, ncx, ncz,
                                                          slow.data(),
                                                          sweepSlowness.data()); 
        } 
        else if (is == 2)
        {
            slownessToSweepSlowness<double, SweepNumber2D::SWEEP3>(nLevels, 
                                                          nx, nz, ncx, ncz,
                                                          slow.data(),
                                                          sweepSlowness.data());
        }
        else //if (is == 3)
        {
            slownessToSweepSlowness<double, SweepNumber2D::SWEEP4>(nLevels, 
                                                          nx, nz, ncx, ncz,
                                                          slow.data(),
                                                          sweepSlowness.data());
        }
        int i0, i1, ix, iz;
        int iCell0X, iCell1X, iCell3X = 0;
        int iCell0Z, iCell1Z, iCell2Z, iCell3Z = 0;
        for (int level = 0; level < nLevels; ++level)
        {
            getLevelStartStopIndices(nx, nz, level, &i0, &i1);
            for (int indx = i0; indx < i1; ++indx)
            {
                // Get the neighbors surrounding the grid point
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP1>(
                                          level, indx, nx, nz, &ix, &iz);
                    iCell0X = std::max(0, ix - 1);
                    iCell1X = std::min(ncx - 1, iCell0X + 1);
                    if (ix == 0){iCell1X = 0;}
                    //iCell2X = iCell1X;
                    iCell3X = iCell0X;
                    iCell0Z = std::max(0, iz - 1);
                    iCell1Z = iCell0Z;
                    iCell2Z = sycl::min(ncz - 1, iCell1Z + 1);
                    if (iz == 0){iCell2Z = 0;}
                    iCell3Z = iCell2Z;
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP2>(
                                          level, indx, nx, nz, &ix, &iz);
                    iCell0X = ix;
                    iCell1X = sycl::max(0, iCell0X - 1);
                    if (ix == ncx)
                    {
                        iCell0X = ix - 1;
                        iCell1X = ix - 1;
                    }
                    //iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = sycl::max(0, iz - 1);
                    iCell2Z = sycl::min(ncz - 1, iCell0Z + 1);
                    if (iz == 0){iCell2Z = 0;}
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP3>(
                                          level, indx, nx, nz, &ix, &iz);
                    iCell0X = sycl::max(0, ix - 1);
                    iCell1X = sycl::min(ncx - 1, iCell0X + 1);
                    if (ix == 0){iCell1X = 0;}
                    //iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = iz;
                    iCell2Z = sycl::max(0, iCell0Z - 1);
                    if (iz == ncz)
                    {
                        iCell0Z = iz - 1;
                        iCell2Z = iz - 1;
                    }
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }
                else // sweep == SweepNumber2D::SWEEP4 
                {
                    sweepLevelIndexToGrid<SweepNumber2D::SWEEP4>(
                                          level, indx, nx, nz, &ix, &iz);
                    iCell0X = ix;
                    iCell1X = sycl::max(0, iCell0X - 1);
                    if (ix == ncx)
                    {
                        iCell0X = ix - 1;
                        iCell1X = ix - 1;
                    }
                    //iCell2X = iCell1X;
                    iCell3X = iCell0X;

                    iCell0Z = iz;
                    iCell2Z = sycl::max(0, iCell0Z - 1);
                    if (iz == ncz)
                    {
                        iCell0Z = iz - 1;
                        iCell2Z = iz - 1;
                    }
                    iCell1Z = iCell0Z;
                    iCell3Z = iCell2Z;
                }

                auto iCell0 = gridToIndex(ncx, iCell0X, iCell0Z);
                auto iCell1 = gridToIndex(ncx, iCell1X, iCell1Z);
                //auto iCell2 = gridToIndex(ncx, iCell2X, iCell2Z);
                auto iCell3 = gridToIndex(ncx, iCell3X, iCell3Z);
                //std::cout << slow[iCell1] << " " << slow[iCell2] << " " << slow[iCell3] << " " << slow[iCell4] << std::endl; 
                int iDst = indx  - i0; 
                EXPECT_NEAR(slow.at(iCell0),
                            sweepSlowness[level].s0[iDst], 1.e-10);
                EXPECT_NEAR(slow.at(iCell1),
                            sweepSlowness[level].s1[iDst], 1.e-10);
                //EXPECT_NEAR(slow.at(iCell2),
                //            sweepSlowness[level].s2[iDst], 1.e-10);
                EXPECT_NEAR(slow.at(iCell3),
                            sweepSlowness[level].s3[iDst], 1.e-10);
            } // Loop on points in level
        }
    }
}

//----------------------------------------------------------------------------//

TEST(Solver2D, loopLimits)
{
    int nx = 10;
    int nz = 11;
    std::vector<int> tRef(4*nx*nz, 0);
    int offset = 0;
    for (int iz = 1; iz < nz; ++iz)
    {
        for (int ix = 1; ix < nx; ++ix)
        {
            tRef.at(offset + gridToIndex(nx, ix, iz)) = 1;
        }
    }

    offset = nx*nz;
    for (int iz = 1; iz < nz; ++iz)
    {
        for (int ix = nx-2; ix >= 0; --ix)
        {
            tRef.at(offset + gridToIndex(nx, ix, iz)) = 2;
        }
    }

    offset = 2*nx*nz;
    for (int iz = nz - 2; iz >= 0; --iz)
    {
        for (int ix = 1; ix < nx; ++ix)
        {
            tRef.at(offset + gridToIndex(nx, ix, iz)) = 3;
        }
    }

    offset = 3*nx*nz;
    for (int iz = nz - 2; iz >= 0; --iz)
    {
        for (int ix = nx - 2; ix >= 0; --ix)
        { 
            tRef.at(offset + gridToIndex(nx, ix, iz)) = 4;
        }
    }

    std::vector<int> t(tRef.size(), 0);
    int ix0, ix1, iz0, iz1, ixDir, izDir;
    for (int is = 0; is < 4; ++is)
    {
        offset = is*nx*nz;
        if (is == 0)
        {
            getLoopLimits<SweepNumber2D::SWEEP1>(nx, nz,
                                                 &ix0, &iz0,
                                                 &ix1, &iz1,
                                                 &ixDir, &izDir);
        }
        else if (is == 1)
        {
            getLoopLimits<SweepNumber2D::SWEEP2>(nx, nz,
                                                 &ix0, &iz0,
                                                 &ix1, &iz1,
                                                 &ixDir, &izDir);
        }
        else if (is == 2)
        {
            getLoopLimits<SweepNumber2D::SWEEP3>(nx, nz,
                                                 &ix0, &iz0,
                                                 &ix1, &iz1,
                                                 &ixDir, &izDir);
        }
        else // if (is == 3)
        {
            getLoopLimits<SweepNumber2D::SWEEP4>(nx, nz,
                                                 &ix0, &iz0,
                                                 &ix1, &iz1,
                                                 &ixDir, &izDir);
        }
        for (int iz = iz0; iz != iz1; iz = iz + izDir)
        {
            for (int ix = ix0; ix != ix1; ix = ix + ixDir)
            {
                t.at(offset + gridToIndex(nx, ix, iz)) = is + 1;
            }
        }
    }
    // Ensure all nodes are visited
    int dTmax = 0;
    for (int i = 0; i < static_cast<int> (t.size()); ++i)
    {
        dTmax = std::max(dTmax, std::abs(t[i] - tRef[i]));
    }
    EXPECT_EQ(dTmax, 0);
}

TEST(Solver2D, maskFSMNodes)
{
    int nx = 11;
    int nz = 10;
    std::vector<int8_t> refUpdateNodes(nx*nz), updateNodes(nx*nz);
    for (int is = 0; is < 4; ++is)
    {
        auto sweep = static_cast<SweepNumber2D> (is);
        std::fill(refUpdateNodes.begin(), refUpdateNodes.end(), UPDATE_NODE);
        for (int iz = 0; iz < nz; ++iz)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                auto indx = gridToIndex(nx, ix, iz);
                if (sweep == SweepNumber2D::SWEEP1)
                {
                    if (ix == 0 || iz == 0)
                    {
                        refUpdateNodes.at(indx) = BOUNDARY_NODE;
                    }
                }
                else if (sweep == SweepNumber2D::SWEEP2)
                {
                    if (ix == nx - 2 || iz == 0)
                    {
                        refUpdateNodes.at(indx) = BOUNDARY_NODE;
                    }
                }
                else if (sweep == SweepNumber2D::SWEEP3)
                {
                    if (ix == 0 || iz == nz - 2)
                    {
                        refUpdateNodes.at(indx) = BOUNDARY_NODE;
                    }
                }
                else // (sweep == SweepNumber2D::SWEEP4)
                {
                    if (ix == nx - 2 || iz == nz - 2)
                    {
                        refUpdateNodes.at(indx) = BOUNDARY_NODE;
                    }
                }
            } // Loop on x
        } // Loop on z
        if (is == 0)
        {
            setPreliminaryUpdateNodes<SweepNumber2D::SWEEP1>(
                nx, nz, updateNodes.data());
        }
        else if (is == 1)
        {
            setPreliminaryUpdateNodes<SweepNumber2D::SWEEP2>(
                nx, nz, updateNodes.data());
        }
        else if (is == 2)
        {
            setPreliminaryUpdateNodes<SweepNumber2D::SWEEP3>(
                nx, nz, updateNodes.data());
        }
        else if (is == 3)
        {
            setPreliminaryUpdateNodes<SweepNumber2D::SWEEP4>(
                nx, nz, updateNodes.data());
        }
        int dmax = 0;
        for (int i = 0; i < static_cast<int> (updateNodes.size()); ++i)
        {
            if (refUpdateNodes[i] != updateNodes[i]){dmax = 1;}
        }
        EXPECT_EQ(dmax, 0);
    } // Loop on sweeps
}

//----------------------------------------------------------------------------//

TEST(Solver2D, solveHomogeneous)
{
    int nx = 21;
    int nz = 20;
    double dx = 100;
    double dz = 100.05;
    double vel = 5000;
    double x0 = 0;
    double z0 = 0;
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double zSrc = z0 + dz*(nz/2);
    int nSweeps = 0;
    int nEps = 3;
    auto solverAlgorithm = EikonalXX::SolverAlgorithm::FAST_SWEEPING_METHOD;

    // Initialize the geometry
    Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);

    // Initialize the velocity model
    std::vector<double> vConstant(nx*nz, vel);
    Model2D<double> vModel;
    vModel.initialize(geometry);
    vModel.setNodalVelocities(vConstant.size(), vConstant.data(),
                              EikonalXX::Ordering2D::NATURAL);

    // Set the source
    Source2D source;
    source.setGeometry(geometry);
    source.setLocationInX(xSrc);
    source.setLocationInZ(zSrc); 

    // Set the solver options
    SolverOptions options;
    options.setNumberOfSweeps(nSweeps);
    options.setFactoredEikonalEquationSolverRadius(nEps);
    options.setVerbosity(Verbosity::DEBUG);
    options.setAlgorithm(solverAlgorithm);
    // Initialize 
    Solver2D<double> solver;
    EXPECT_NO_THROW(solver.initialize(options, geometry));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_NO_THROW(solver.setVelocityModel(vModel));
    EXPECT_NO_THROW(solver.setSource(std::pair(xSrc, zSrc)));
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_TRUE(solver.haveSource());
    // Solve
//    EXPECT_NO_THROW(solver.solve());
/*
    EXPECT_TRUE(solver.haveTravelTimeField());
    auto tEst = solver.getTravelTimeField();
    std::ofstream ofl("ttimes.txt");
    for (int iz=0; iz<nz; ++iz)
    {
        for (int ix=0; ix<nx; ++ix)
        {
            ofl << x0 + ix*dx << " " << z0 + iz*dz << " "
                << tEst[gridToIndex(nx, ix, iz)] << std::endl;
        }
        ofl << std::endl;
    }
    ofl.close();
*/
} 

TEST(Solver2D, Increment)
{
    int nx = 21;
    int nz = 20;
    double dx = 100; 
    double dz = 100.05;
    double x0 = 1; 
    double z0 = 2; 
    double xSrc = x0 + dx*(nx/2.) + dx/4.;
    double zSrc = z0 + dz*(nz/2.);
    int nSweeps = 0; 
    int nEps = 3; 
    auto solverAlgorithm = EikonalXX::SolverAlgorithm::FAST_SWEEPING_METHOD;

    // Initialize the geometry
    Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);

    // Initialize the velocity model
    std::vector<double> vIncrement(nx*nz, 0);
    int ic = 0;
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            ic = ic + 1;
            auto indx = gridToIndex(nx, ix, iz);
            vIncrement.at(indx) = ic;//10*(indx + 1) + 1;
        }
    }
    //getchar();
    Model2D<double> vModel;
    vModel.initialize(geometry);
    vModel.setNodalVelocities(vIncrement.size(), vIncrement.data(),
                              EikonalXX::Ordering2D::NATURAL);

    // Set the source
    Source2D source;
    source.setGeometry(geometry);
    source.setLocationInX(xSrc);
    source.setLocationInZ(zSrc);

    // Set the solver options
    SolverOptions options;
    options.setNumberOfSweeps(nSweeps);
    options.setFactoredEikonalEquationSolverRadius(nEps);
    options.setVerbosity(Verbosity::DEBUG);
    options.setAlgorithm(solverAlgorithm);
    // Initialize 
    Solver2D<double> solver;
    EXPECT_NO_THROW(solver.initialize(options, geometry));
    EXPECT_TRUE(solver.isInitialized());
    EXPECT_NO_THROW(solver.setVelocityModel(vModel));
    EXPECT_NO_THROW(solver.setSource(std::pair(xSrc, zSrc)));
    EXPECT_TRUE(solver.haveVelocityModel());
    EXPECT_TRUE(solver.haveSource());
    // Solve
    EXPECT_NO_THROW(solver.solve());
}

/*
TEST(Solver2D, finiteDifferenceShift)
{
    int shiftTx, shiftTz, shiftVx, shiftVz;
    // +x and +z
    getSweepFiniteDifferenceShifts(1, &shiftTx, &shiftTz, &shiftVx, &shiftVz);
    EXPECT_EQ(shiftTx, 1);
    EXPECT_EQ(shiftTz, 1);
    EXPECT_EQ(shiftVx, 1);
    EXPECT_EQ(shiftVz, 1);
    // -x and +z
    getSweepFiniteDifferenceShifts(2, &shiftTx, &shiftTz, &shiftVx, &shiftVz);
    EXPECT_EQ(shiftTx,-1);
    EXPECT_EQ(shiftTz, 1);
    EXPECT_EQ(shiftVx, 0);
    EXPECT_EQ(shiftVz, 1);
    // +x and -z
    getSweepFiniteDifferenceShifts(3, &shiftTx, &shiftTz, &shiftVx, &shiftVz);
    EXPECT_EQ(shiftTx, 1);
    EXPECT_EQ(shiftTz,-1);
    EXPECT_EQ(shiftVx, 1);
    EXPECT_EQ(shiftVz, 0);
    // -x and -z
    getSweepFiniteDifferenceShifts(4, &shiftTx, &shiftTz, &shiftVx, &shiftVz);
    EXPECT_EQ(shiftTx,-1);
    EXPECT_EQ(shiftTz,-1);
    EXPECT_EQ(shiftVx, 0);
    EXPECT_EQ(shiftVz, 0);

}
*/

}
