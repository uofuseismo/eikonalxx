#include "eikonalxx/graph3d.hpp"
#include "private/grid.hpp"
#include <gtest/gtest.h>
namespace
{
using namespace EikonalXX;

TEST(Graph, indexToGrid)
{
    int nx = 44;
    int ny = 53;
    int nz = 23;
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                auto indx = gridToIndex(nx, ny, ix, iy, iz);
                EXPECT_EQ(indx, iz*nx*ny + iy*nx + ix);
                int jx, jy, jz;
                indexToGrid(indx, nx, ny, &jx, &jy, &jz); 
                EXPECT_EQ(ix, jx);
                EXPECT_EQ(iy, jy);
                EXPECT_EQ(iz, jz);
            }
        }
    }
}

TEST(Graph, numberOfLevels3D)
{
    int nx = 5;  
    int ny = 3;
    int nz = 7;  
    int nLevels = nx + ny + nz - 2;  
    EXPECT_EQ(computeNumberOfLevels(nx, ny, nz), nLevels);
}

TEST(Graph, Graph3D)
{
    int nx = 4;
    int ny = 3;
    int nz = 5;
    Graph3D<SweepNumber3D::SWEEP1> graph;
    EXPECT_NO_THROW(graph.initialize(nx, ny, nz));
    EXPECT_TRUE(graph.isInitialized());

}
}
