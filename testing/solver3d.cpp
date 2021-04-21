#include <iostream>
#include <fstream>
#include <vector>
#include "eikonalxx/solver3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "eikonalxx/source3d.hpp"
#include "eikonalxx/model3d.hpp"
//#include "private/solverUtilities2d.hpp"
#include "private/solverUtilities3d.hpp"
//#include "private/solver3d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(Solver3D, analyticSolutions)
{
    double dx = 100.;
    double dy = 150;
    double dz = 200.;
    double xSrc = 500;
    double ySrc = 750;
    double zSrc = 1000; 
    double slowness = 1./5000.;
    double t, dtdx, dtdy, dtdz;
    
    int ix = 1;
    int iy = 1;
    int iz = 1;
    auto x = ix*dx;
    auto y = iy*dy;
    auto z = iz*dz;
    double distance = std::sqrt((x - xSrc)*(x - xSrc) 
                              + (y - ySrc)*(y - ySrc)
                              + (z - zSrc)*(z - zSrc));
    double tRef = distance*slowness;
    // dT/dx, dT/dy, and dT/dz are really just apparent slownesses
    // so compute normal vector in that direction then scale by slowness.
    auto nx = (x - xSrc)/distance;
    auto ny = (y - ySrc)/distance;
    auto nz = (z - zSrc)/distance;
    auto sx = nx*slowness;
    auto sy = ny*slowness;
    auto sz = nz*slowness;

    t = computeAnalyticalTravelTime(ix, iy, iz, dx, dy, dz,
                                    xSrc, ySrc, zSrc, slowness);
    EXPECT_NEAR(t, tRef, 1.e-14);

    computeAnalyticalTravelTime(ix, iy, iz,
                                dx, dy, dz,
                                xSrc, ySrc, zSrc,
                                slowness,
                                &t, &dtdx, &dtdy, &dtdz);
    EXPECT_NEAR(t,  tRef, 1.e-14);
    EXPECT_NEAR(sx, dtdx, 1.e-14);
    EXPECT_NEAR(sy, dtdy, 1.e-14);
    EXPECT_NEAR(sz, dtdz, 1.e-14); 
}

TEST(Solver3D, sweepSigns)
{
    int ixShift, iyShift, izShift;
    int signX, signY, signZ;
    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP1>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, -1); EXPECT_EQ(signX,  1);
    EXPECT_EQ(iyShift, -1); EXPECT_EQ(signY,  1);
    EXPECT_EQ(izShift, -1); EXPECT_EQ(signZ,  1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP2>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, +1); EXPECT_EQ(signX, -1);
    EXPECT_EQ(iyShift, -1); EXPECT_EQ(signY,  1);
    EXPECT_EQ(izShift, -1); EXPECT_EQ(signZ,  1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP3>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, -1); EXPECT_EQ(signX,  1);
    EXPECT_EQ(iyShift, +1); EXPECT_EQ(signY, -1);
    EXPECT_EQ(izShift, -1); EXPECT_EQ(signZ,  1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP4>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, +1); EXPECT_EQ(signX, -1); 
    EXPECT_EQ(iyShift, +1); EXPECT_EQ(signY, -1); 
    EXPECT_EQ(izShift, -1); EXPECT_EQ(signZ,  1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP5>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, -1); EXPECT_EQ(signX,  1);
    EXPECT_EQ(iyShift, -1); EXPECT_EQ(signY,  1);
    EXPECT_EQ(izShift, +1); EXPECT_EQ(signZ, -1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP6>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, +1); EXPECT_EQ(signX, -1);
    EXPECT_EQ(iyShift, -1); EXPECT_EQ(signY,  1);
    EXPECT_EQ(izShift, +1); EXPECT_EQ(signZ, -1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP7>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, -1); EXPECT_EQ(signX,  1);
    EXPECT_EQ(iyShift, +1); EXPECT_EQ(signY, -1);
    EXPECT_EQ(izShift, +1); EXPECT_EQ(signZ, -1);

    getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP8>(
        &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
    EXPECT_EQ(ixShift, +1); EXPECT_EQ(signX, -1); 
    EXPECT_EQ(iyShift, +1); EXPECT_EQ(signY, -1);
    EXPECT_EQ(izShift, +1); EXPECT_EQ(signZ, -1);
}

TEST(Solver3D, sweepToSlownessIndices)
{
    int nx = 25;
    int ny = 26;
    int nz = 27;
    auto nCellX = nx - 1;
    auto nCellY = ny - 1;
    auto nCellZ = nz - 1;
    auto nCell = nCellX*nCellY*nCellZ;
    int ixShift, iyShift, izShift;
    int i0, i1, i2, i3, i4, i5, i7;
    int signX, signY, signZ;
    int i0Ref, i1Ref, i2Ref, i3Ref, i4Ref, i5Ref, i7Ref;
    for (int is = 0; is < 8; ++is)
    {
        int sgnvx = 1;
        int sgnvy = 1;
        int sgnvz = 1;
        if (is == 0)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP1>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 1)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP2>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvx = 0;
        }
        else if (is == 2)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP3>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvy = 0;
        }
        else if (is == 3)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP4>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvx = 0;
            sgnvy = 0;
        }
        else if (is == 4)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP5>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvz = 0;
        }
        else if (is == 5)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP6>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvx = 0;
            sgnvz = 0;
        }
        else if (is == 6)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP7>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvy = 0;
            sgnvz = 0;
        }
        else
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP8>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
            sgnvx = 0;
            sgnvy = 0;
            sgnvz = 0;
        }
        // Loop on grid
        for (int iz = 0; iz < nz; ++iz)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    auto iCellX = std::min(nCellX - 1, std::max(0, ix - sgnvx));
                    auto iCellY = std::min(nCellY - 1, std::max(0, iy - sgnvy));
                    auto iCellZ = std::min(nCellZ - 1, std::max(0, iz - sgnvz));
                    i0Ref = gridToIndex(nCellX, nCellY,
                                        iCellX, iCellY, iCellZ);
                    i1Ref = gridToIndex(nCellX, nCellY,
                            std::max(0, std::min(nCellX - 1, iCellX - ixShift)),
                            iCellY,
                            iCellZ);
                    i2Ref = gridToIndex(nCellX, nCellY,
                            std::max(0, std::min(nCellX - 1, iCellX - ixShift)),
                            std::max(0, std::min(nCellY - 1, iCellY - iyShift)),
                            iCellZ);
                    i3Ref = gridToIndex(nCellX, nCellY,
                            iCellX,
                            std::max(0, std::min(nCellY - 1, iCellY - iyShift)), 
                            iCellZ);
                    i4Ref = gridToIndex(nCellX, nCellY,
                            iCellX,
                            iCellY,
                            std::max(0, std::min(nCellZ - 1, iCellZ - izShift)));
                    i5Ref = gridToIndex(nCellX, nCellY,
                            std::max(0, std::min(nCellX - 1, iCellX - ixShift)),
                            iCellY,
                            std::max(0, std::min(nCellZ - 1, iCellZ - izShift)));
                    i7Ref = gridToIndex(nCellX, nCellY,
                            iCellX,
                            std::max(0, std::min(nCellY - 1, iCellY - iyShift)),
                            std::max(0, std::min(nCellZ - 1, iCellZ - izShift)));
                    if (is == 0)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP1>(
                                                 ix, iy, iz, 
                                                 nCellX, nCellY, nCellZ, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 1)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP2>(
                                                 ix, iy, iz, 
                                                 nCellX, nCellY, nCellZ, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 2)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP3>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 3)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP4>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 4)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP5>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 5)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP6>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else if (is == 6)
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP7>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    else
                    {
                        gridToSurroundingSlownessIndices<SweepNumber3D::SWEEP8>(
                                                 ix, iy, iz,
                                                 nCellX, nCellY, nCellZ,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i7);
                    }
                    EXPECT_TRUE(i0 >= 0 && i0 < nCell);
                    EXPECT_TRUE(i1 >= 0 && i1 < nCell);
                    EXPECT_TRUE(i2 >= 0 && i2 < nCell);
                    EXPECT_TRUE(i3 >= 0 && i3 < nCell);
                    EXPECT_TRUE(i4 >= 0 && i4 < nCell);
                    EXPECT_TRUE(i5 >= 0 && i5 < nCell);
                    //EXPECT_TRUE(i6 >= 0 && i6 < nCell);
                    EXPECT_TRUE(i7 >= 0 && i7 < nCell);
                    EXPECT_EQ(i0, i0Ref);
                    EXPECT_EQ(i1, i1Ref);
                    EXPECT_EQ(i2, i2Ref);
                    EXPECT_EQ(i3, i3Ref);
                    EXPECT_EQ(i4, i4Ref);
                    EXPECT_EQ(i5, i5Ref);
                    EXPECT_EQ(i7, i7Ref);
                }
            }
        }
    }
}

TEST(Solver3D, sweepToTravelTimeIndices)
{
    int nx = 32;
    int ny = 43;
    int nz = 56;
    auto nGrid = nx*ny*nz;
    int ixShift, iyShift, izShift;
    int signX, signY, signZ;
    int i0, i1, i2, i3, i4, i5, i6, i7;
    int i0Ref, i1Ref, i2Ref, i3Ref, i4Ref, i5Ref, i6Ref, i7Ref;
    for (int is = 0; is < 8; ++is)
    {
        if (is == 0)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP1>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 1)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP2>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 2)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP3>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 3)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP4>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 4)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP5>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 5)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP6>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else if (is == 6)
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP7>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        else
        {
            getSweepFiniteDifferenceSigns<SweepNumber3D::SWEEP8>(
                &ixShift, &iyShift, &izShift, &signX, &signY, &signZ);
        }
        // Loop on grid
        for (int iz = 0; iz < nz; ++iz)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    i0Ref = gridToIndex(nx, ny, ix, iy, iz);
                    i1Ref = gridToIndex(nx, ny, 
                                   std::max(0, std::min(nx - 1, ix + ixShift)),
                                   iy,
                                   iz); 
                    i2Ref = gridToIndex(nx, ny,
                                   std::max(0, std::min(nx - 1, ix + ixShift)),
                                   std::max(0, std::min(ny - 1, iy + iyShift)),
                                   iz);
                    i3Ref = gridToIndex(nx, ny,
                                   ix,
                                   std::max(0, std::min(ny - 1, iy + iyShift)),
                                   iz);
                    i4Ref = gridToIndex(nx ,ny,
                                   ix,
                                   iy,
                                   std::max(0, std::min(nz - 1, iz + izShift)));
                    i5Ref = gridToIndex(nx, ny, 
                                   std::max(0, std::min(nx - 1, ix + ixShift)),
                                   iy, 
                                   std::max(0, std::min(nz - 1, iz + izShift)));
                    i6Ref = gridToIndex(nx, ny, 
                                   std::max(0, std::min(nx - 1, ix + ixShift)),
                                   std::max(0, std::min(ny - 1, iy + iyShift)),
                                   std::max(0, std::min(nz - 1, iz + izShift)));
                    i7Ref = gridToIndex(nx, ny, 
                                   ix, 
                                   std::max(0, std::min(ny - 1, iy + iyShift)),
                                   std::max(0, std::min(nz - 1, iz + izShift)));
                    if (is == 0)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP1>(
                                                 ix, iy, iz,
                                                 nx, ny, nz,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 1)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP2>(
                                                 ix, iy, iz, 
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 2)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP3>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 3)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP4>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 4)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP5>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 5)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP6>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 6)
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP7>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else
                    {
                        gridToSurroundingTravelTimeIndices<SweepNumber3D::SWEEP8>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    EXPECT_TRUE(i0 >= 0 && i0 < nGrid);
                    EXPECT_TRUE(i1 >= 0 && i1 < nGrid);
                    EXPECT_TRUE(i2 >= 0 && i2 < nGrid);
                    EXPECT_TRUE(i3 >= 0 && i3 < nGrid);
                    EXPECT_TRUE(i4 >= 0 && i4 < nGrid);
                    EXPECT_TRUE(i5 >= 0 && i5 < nGrid);
                    EXPECT_TRUE(i6 >= 0 && i6 < nGrid);
                    EXPECT_TRUE(i7 >= 0 && i7 < nGrid);
                    EXPECT_EQ(i0Ref, i0);
                    EXPECT_EQ(i1Ref, i1);
                    EXPECT_EQ(i2Ref, i2);
                    EXPECT_EQ(i3Ref, i3);
                    EXPECT_EQ(i4Ref, i4);
                    EXPECT_EQ(i5Ref, i5);
                    EXPECT_EQ(i6Ref, i6);
                    EXPECT_EQ(i7Ref, i7);
                    //std::cout << "(" << ix << "," << iy << "," << iz << ") "
                    //          << i0 << " " << i1 << " " << i2 << " " << i3 << " "
                    //          << i4 << " " << i5 << " " << i6 << " " << i7
                    //          <<  std::endl;
                }
            }
        }
    }
}

TEST(Solver3D, loopLimits)
{
    int nGridX = 10;
    int nGridY = 11;
    int nGridZ = 9;
    int ix0, iy0, iz0, ix1, iy1, iz1, ixDir, iyDir, izDir;
    getLoopLimits<SweepNumber3D::SWEEP1>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, 1); EXPECT_EQ(ix1, nGridX); EXPECT_EQ(ixDir, 1);
    EXPECT_EQ(iy0, 1); EXPECT_EQ(iy1, nGridY); EXPECT_EQ(iyDir, 1);
    EXPECT_EQ(iz0, 1); EXPECT_EQ(iz1, nGridZ); EXPECT_EQ(izDir, 1);

    getLoopLimits<SweepNumber3D::SWEEP2>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, nGridX - 2); EXPECT_EQ(ix1,-1);      EXPECT_EQ(ixDir,-1);
    EXPECT_EQ(iy0, 1);          EXPECT_EQ(iy1, nGridY); EXPECT_EQ(iyDir, 1); 
    EXPECT_EQ(iz0, 1);          EXPECT_EQ(iz1, nGridZ); EXPECT_EQ(izDir, 1);

    getLoopLimits<SweepNumber3D::SWEEP3>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, 1);          EXPECT_EQ(ix1, nGridX); EXPECT_EQ(ixDir, 1); 
    EXPECT_EQ(iy0, nGridY - 2); EXPECT_EQ(iy1,-1);      EXPECT_EQ(iyDir,-1);
    EXPECT_EQ(iz0, 1);          EXPECT_EQ(iz1, nGridZ); EXPECT_EQ(izDir, 1); 

    getLoopLimits<SweepNumber3D::SWEEP4>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, nGridX - 2); EXPECT_EQ(ix1,-1);      EXPECT_EQ(ixDir,-1);
    EXPECT_EQ(iy0, nGridY - 2); EXPECT_EQ(iy1,-1);      EXPECT_EQ(iyDir,-1);
    EXPECT_EQ(iz0, 1);          EXPECT_EQ(iz1, nGridZ); EXPECT_EQ(izDir, 1);

    getLoopLimits<SweepNumber3D::SWEEP5>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, 1);          EXPECT_EQ(ix1, nGridX); EXPECT_EQ(ixDir, 1); 
    EXPECT_EQ(iy0, 1);          EXPECT_EQ(iy1, nGridY); EXPECT_EQ(iyDir, 1); 
    EXPECT_EQ(iz0, nGridZ - 2); EXPECT_EQ(iz1,-1);      EXPECT_EQ(izDir,-1); 

    getLoopLimits<SweepNumber3D::SWEEP6>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, nGridX - 2); EXPECT_EQ(ix1,-1);      EXPECT_EQ(ixDir,-1);
    EXPECT_EQ(iy0, 1);          EXPECT_EQ(iy1, nGridY); EXPECT_EQ(iyDir, 1); 
    EXPECT_EQ(iz0, nGridZ - 2); EXPECT_EQ(iz1,-1);      EXPECT_EQ(izDir,-1);

    getLoopLimits<SweepNumber3D::SWEEP7>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, 1);          EXPECT_EQ(ix1, nGridX); EXPECT_EQ(ixDir, 1); 
    EXPECT_EQ(iy0, nGridY - 2); EXPECT_EQ(iy1,-1);      EXPECT_EQ(iyDir,-1);
    EXPECT_EQ(iz0, nGridZ - 2); EXPECT_EQ(iz1,-1);      EXPECT_EQ(izDir,-1);

    getLoopLimits<SweepNumber3D::SWEEP8>(nGridX, nGridY, nGridZ,
                                         &ix0, &iy0, &iz0,
                                         &ix1, &iy1, &iz1,
                                         &ixDir, &iyDir, &izDir);
    EXPECT_EQ(ix0, nGridX - 2); EXPECT_EQ(ix1,-1);      EXPECT_EQ(ixDir,-1);
    EXPECT_EQ(iy0, nGridY - 2); EXPECT_EQ(iy1,-1);      EXPECT_EQ(iyDir,-1);
    EXPECT_EQ(iz0, nGridZ - 2); EXPECT_EQ(iz1,-1);      EXPECT_EQ(izDir,-1);
}

TEST(Solver3D, permuteGrid)
{
    int nx = 4;
    int ny = 5; 
    int nz = 6;
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                int jx, jy, jz;
                permuteGrid<SweepNumber3D::SWEEP1>(ix, iy, iz, nx, ny, nz,
                                                   &jx, &jy, &jz);
                EXPECT_EQ(ix, jx);
                EXPECT_EQ(iy, jy);
                EXPECT_EQ(iz, jz);
                permuteGrid<SweepNumber3D::SWEEP2>(ix, iy, iz, nx, ny, nz, 
                                                   &jx, &jy, &jz);
                EXPECT_EQ(nx - 1 - ix, jx);
                EXPECT_EQ(iy, jy);
                EXPECT_EQ(iz, jz);
                permuteGrid<SweepNumber3D::SWEEP3>(ix, iy, iz, nx, ny, nz, 
                                                   &jx, &jy, &jz);
                EXPECT_EQ(ix, jx);
                EXPECT_EQ(ny - 1 - iy, jy);
                EXPECT_EQ(iz, jz);
                permuteGrid<SweepNumber3D::SWEEP4>(ix, iy, iz, nx, ny, nz, 
                                                   &jx, &jy, &jz);
                EXPECT_EQ(nx - 1 - ix, jx); 
                EXPECT_EQ(ny - 1 - iy, jy);
                EXPECT_EQ(iz, jz);

                permuteGrid<SweepNumber3D::SWEEP5>(ix, iy, iz, nx, ny, nz, 
                                                   &jx, &jy, &jz);
                EXPECT_EQ(ix, jx);
                EXPECT_EQ(iy, jy);
                EXPECT_EQ(nz - 1 - iz, jz);
                permuteGrid<SweepNumber3D::SWEEP6>(ix, iy, iz, nx, ny, nz,
                                                   &jx, &jy, &jz);
                EXPECT_EQ(nx - 1 - ix, jx);
                EXPECT_EQ(iy, jy);
                EXPECT_EQ(nz - 1 - iz, jz);
                permuteGrid<SweepNumber3D::SWEEP7>(ix, iy, iz, nx, ny, nz,
                                                   &jx, &jy, &jz);
                EXPECT_EQ(ix, jx);
                EXPECT_EQ(ny - 1 - iy, jy);
                EXPECT_EQ(nz - 1 - iz, jz);
                permuteGrid<SweepNumber3D::SWEEP8>(ix, iy, iz, nx, ny, nz,
                                                   &jx, &jy, &jz);
                EXPECT_EQ(nx - 1 - ix, jx);
                EXPECT_EQ(ny - 1 - iy, jy);
                EXPECT_EQ(nz - 1 - iz, jz);
            }
        }
    }
}

TEST(Solver3D, isSweepBoundaryNode)
{
    int nx = 5;
    int ny = 6;
    int nz = 7;
    for (int is = 0; is < 8; ++is)
    {
        for (int iz =0 ; iz < nz; ++iz)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    bool lb = false;
                    bool lbRef = false;
                    if (is == 0)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP1>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == 0 || iy == 0 || iz == 0){lbRef = true;} 
                    }
                    else if (is == 1)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP2>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == nx - 1 || iy == 0 || iz == 0){lbRef = true;}
                    }
                    else if (is == 2)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP3>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == 0 || iy == ny - 1 || iz == 0){lbRef = true;}
                    }
                    else if (is == 3)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP4>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == nx - 1 || iy == ny - 1 || iz == 0)
                        {
                            lbRef = true;
                        }
                    }
                    else if (is == 4)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP5>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == 0 || iy == 0 || iz == nz - 1){lbRef = true;} 
                    }
                    else if (is == 5)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP6>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == nx - 1 || iy == 0 || iz == nz - 1)
                        {
                            lbRef = true;
                        }
                    }
                    else if (is == 6)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP7>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == 0 || iy == ny - 1 || iz == nz - 1)
                        {
                            lbRef = true;
                        }
                    }
                    else if (is == 7)
                    {
                        lb = isSweepBoundaryNode<SweepNumber3D::SWEEP8>(
                                 ix, iy, iz, nx, ny, nz);
                        if (ix == nx - 1 || iy == ny - 1 || iz == nz - 1)
                        {
                            lbRef = true;
                        }
                    }
                    EXPECT_TRUE(lb == lbRef);
                }
            }
        }
    }
}

TEST(Solver3D, solve3d)
{
std::cout.precision(10);
    int nx = 9;
    int ny = 10;
    int nz = 11;
    double dx = 100;
    double dy = dx;
    double dz = dx;
    double x0 = 5;
    double y0 = 6;
    double z0 = 7;
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double ySrc = y0 + dy*(ny/2) - dy/4;
    double zSrc = z0 + dz*(nz/4);
    double vConst = 3000;
    int nSweeps = 1;
    int nEps =-1;
nEps = 500; //500;
    auto solverAlgorithm = EikonalXX::SolverAlgorithm::FAST_SWEEPING_METHOD;

    SolverOptions options;
    options.setVerbosity(Verbosity::DEBUG);
    options.setNumberOfSweeps(nSweeps);
    options.setSphericalSolverRadius(nEps);
    options.setAlgorithm(solverAlgorithm);

    // Set the geometry
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

    // Set the velocity model
    std::vector<double> vConstant(nx*ny*nz, vConst);
    Model3D<double> vModel;
    vModel.initialize(geometry);
    vModel.setNodalVelocities(vConstant.size(), vConstant.data(),
                              EikonalXX::Ordering3D::NATURAL);

    // Set the source
    Source3D source;
    source.setGeometry(geometry);
    source.setLocationInX(xSrc);
    source.setLocationInY(ySrc);
    source.setLocationInZ(zSrc);

    // Initialize the solver
    Solver3D<double> solver;
    EXPECT_NO_THROW(solver.initialize(options, geometry));
    EXPECT_NO_THROW(solver.setVelocityModel(vModel));
    EXPECT_NO_THROW(solver.setSource(source));
    EXPECT_NO_THROW(solver.solve());
} 

}
