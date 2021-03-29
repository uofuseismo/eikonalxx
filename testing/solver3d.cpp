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

TEST(Solver3D, sweepToGridTravelTimeIndices)
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
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP1>(
                                                 ix, iy, iz,
                                                 nx, ny, nz,
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 1)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP2>(
                                                 ix, iy, iz, 
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 2)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP3>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 3)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP4>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 4)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP5>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 5)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP6>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else if (is == 6)
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP7>(
                                                 ix, iy, iz,
                                                 nx, ny, nz, 
                                                 &i0, &i1, &i2, &i3,
                                                 &i4, &i5, &i6, &i7);
                    }
                    else
                    {
                        gridToSurroundingTravelTimes<SweepNumber3D::SWEEP8>(
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

TEST(Solver3D, solve3d)
{
    int nx = 9;
    int ny = 10;
    int nz = 11; 
    double dx = 100;
    double dy = dx;
    double dz = dx;
    double x0 = 0;
    double y0 = 0;
    double z0 = 0;
    double xSrc = x0 + dx*(nx/2) + dx/4;
    double ySrc = y0 + dy*(ny/2) - dy/4;
    double zSrc = z0 + dz*(nz/2);
    SolverOptions options;

    // Set the geometry
    Geometry3D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInY(ny);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInY(dy);
    geometry.setGridSpacingInZ(dz); 

    // Set the source
    Source3D source;
    source.setGeometry(geometry);
    source.setLocationInX(xSrc);
    source.setLocationInY(ySrc);
    source.setLocationInZ(zSrc);

    // Initialize the solver
    Solver3D<double> solver;
    solver.initialize(options, geometry);
    solver.setSource(source);
} 

}
