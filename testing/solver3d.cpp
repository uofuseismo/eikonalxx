#include <iostream>
#include <fstream>
#include <vector>
//#include "eikonalxx/solver3d.hpp"
//#include "eikonalxx/geometry3d.hpp"
//#include "eikonalxx/model3d.hpp"
//#include "private/solverUtilities2d.hpp"
#include "private/solverUtilities3d.hpp"
//#include "private/solver3d.hpp"
#include "eikonalxx/solverOptions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;


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

}
