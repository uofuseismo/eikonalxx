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
    int nLevelsRef = 10;
    std::vector<int> nodeInLevelToGridPointRef1{
                                0,
                                1, 4,12,
                                2, 5, 8,13,16,24,
                                3, 6, 9,14,17,20,25,28,36,
                                7,10,15,18,21,26,29,32,37,40,48,
                               11,19,22,27,30,33,38,41,44,49,52,
                               23,31,34,39,42,45,50,53,56,
                               35,43,46,51,54,57,
                               47,55,58,
                               59};
    std::vector<int> nodeInLevelToGridPointRef2{
                                3,
                                2, 7,15,
                                1, 6,11,14,19,27,
                                0, 5,10,13,18,23,26,31,39,
                                4, 9,12,17,22,25,30,35,38,43,51,
                                8,16,21,24,29,34,37,42,47,50,55,
                               20,28,33,36,41,46,49,54,59,
                               32,40,45,48,53,58,
                               44,52,57,
                               56};
    std::vector<int> nodeInLevelToGridPointRef3{
                                8,
                                4, 9,20,
                                0, 5,10,16,21,32,
                                1, 6,11,12,17,22,28,33,44,
                                2, 7,13,18,23,24,29,34,40,45,56,
                                3,14,19,25,30,35,36,41,46,52,57,
                               15,26,31,37,42,47,48,53,58,
                               27,38,43,49,54,59,
                               39,50,55,
                               51};
    std::vector<int> nodeInLevelToGridPointRef4{
                               11,
                                7,10,23,
                                3, 6, 9,19,22,35,
                                2, 5, 8,15,18,21,31,34,47,
                                1, 4,14,17,20,27,30,33,43,46,59,
                                0,13,16,26,29,32,39,42,45,55,58,
                               12,25,28,38,41,44,51,54,57,
                               24,37,40,50,53,56,
                               36,49,52,
                               48};
    std::vector<int> nodeInLevelToGridPointRef5{ // Needs to be sorted
                               48,
                               36,49,52, //49,52,36,
                               24,37,40,50,53,56, //50,53,56,37,40,24,
                               12,25,28,38,41,44,51,54,57, //51,54,57,38,41,44,25,28,12,
                                0,13,16,26,29,32,39,42,45,55,58, //55,58,39,42,45,26,29,32,13,16, 0,
                                1, 4,14,17,20,27,30,33,43,46,59, //59,43,46,27,30,33,14,17,20, 1, 4,
                                2, 5, 8,15,18,21,31,34,47, //47,31,34,15,18,21, 2, 5, 8,
                                3, 6, 9,19,22,35, //35,19,22, 3, 6, 9,
                                7,10,23, //23, 7,10,
                               11};
    std::vector<int> nodeInLevelToGridPointRef6{
                               51,
                               39,50,55,
                               27,38,43,49,54,59,
                               15,26,31,37,42,47,48,53,58,
                                3,14,19,25,30,35,36,41,46,52,57,
                                2, 7,13,18,23,24,29,34,40,45,56,
                                1, 6,11,12,17,22,28,33,44,
                                0, 5,10,16,21,32,
                                4, 9,20,
                                8};
    std::vector<int> nodeInLevelToGridPointRef7{
                               56,
                               44,52,57,
                               32,40,45,48,53,58,
                               20,28,33,36,41,46,49,54,59,
                                8,16,21,24,29,34,37,42,47,50,55,
                                4, 9,12,17,22,25,30,35,38,43,51,
                                0, 5,10,13,18,23,26,31,39,
                                1, 6,11,14,19,27,
                                2, 7,15,
                                3}; 
    std::vector<int> nodeInLevelToGridPointRef8{
                               59,
                               47,55,58,
                               35,43,46,51,54,57,
                               23,31,34,39,42,45,50,53,56,
                               11,19,22,27,30,33,38,41,44,49,52,
                                7,10,15,18,21,26,29,32,37,40,48,
                                3, 6, 9,14,17,20,25,28,36,
                                2, 5, 8,13,16,24,
                                1, 4,12,
                                0};
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef1.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef2.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef3.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef4.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef5.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef6.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef7.size());
    ASSERT_EQ(static_cast<size_t> (nx*ny*nz),
              nodeInLevelToGridPointRef8.size());

    int maxLevelSizeRef = 11;
    std::vector<int> levelStartPtrRef{0,
                                      1,
                                      4,
                                      10,
                                      19,
                                      30,
                                      41,
                                      50,
                                      56,
                                      59,
                                      60};
                                            
    Graph3D<SweepNumber3D::Sweep1> graph1;
    EXPECT_NO_THROW(graph1.initialize(nx, ny, nz));
    EXPECT_TRUE(graph1.isInitialized());
    EXPECT_EQ(nx*ny*nz, graph1.getNumberOfGridPoints());
    EXPECT_EQ(nLevelsRef, graph1.getNumberOfLevels());
    // Test level pointer
    auto levelStartPtr = graph1.getLevelStartPointer();
    for (int level = 0; level < graph1.getNumberOfLevels() + 1; ++level)
    {
        EXPECT_EQ(levelStartPtrRef[level], levelStartPtr[level]);
        if (level < graph1.getNumberOfLevels())
        {
            EXPECT_EQ(levelStartPtrRef[level + 1] - levelStartPtrRef[level],
                      graph1.getNumberOfNodesInLevel(level));
        }
    }
    EXPECT_EQ(maxLevelSizeRef, graph1.getMaximumLevelSize());
    auto nodeInLevelToGridPointPtr = graph1.getNodeInLevelToIndexPointer();
    auto nodeToXGridPointPtr = graph1.getNodeInLevelToXGridPointPointer();
    auto nodeToYGridPointPtr = graph1.getNodeInLevelToYGridPointPointer();
    auto nodeToZGridPointPtr = graph1.getNodeInLevelToZGridPointPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {
        EXPECT_EQ(nodeInLevelToGridPointRef1[i],
                  nodeInLevelToGridPointPtr[i]);
        int ix, iy, iz;
        indexToGrid(nodeInLevelToGridPointRef1[i], nx, ny, &ix, &iy, &iz);
        EXPECT_EQ(ix, nodeToXGridPointPtr[i]);
        EXPECT_EQ(iy, nodeToYGridPointPtr[i]);
        EXPECT_EQ(iz, nodeToZGridPointPtr[i]);
    } 

    // Sweep direction 2
    Graph3D<SweepNumber3D::Sweep2> graph2;
    EXPECT_NO_THROW(graph2.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph2.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {   
        EXPECT_EQ(nodeInLevelToGridPointRef2[i],
                  nodeInLevelToGridPointPtr[i]);
    }   

    // Sweep direction 3
    Graph3D<SweepNumber3D::Sweep3> graph3;
    EXPECT_NO_THROW(graph3.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph3.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {
        EXPECT_EQ(nodeInLevelToGridPointRef3[i],
                  nodeInLevelToGridPointPtr[i]);
    }   

    // Sweep direction 4
    Graph3D<SweepNumber3D::Sweep4> graph4;
    EXPECT_NO_THROW(graph4.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph4.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {   
        EXPECT_EQ(nodeInLevelToGridPointRef4[i],
                  nodeInLevelToGridPointPtr[i]);
    }

    // Sweep direction 5
    Graph3D<SweepNumber3D::Sweep5> graph5;
    EXPECT_NO_THROW(graph5.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph5.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {
        EXPECT_EQ(nodeInLevelToGridPointRef5[i],
                  nodeInLevelToGridPointPtr[i]);
    }

    // Sweep direction 6
    Graph3D<SweepNumber3D::Sweep6> graph6;
    EXPECT_NO_THROW(graph6.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph6.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {
        EXPECT_EQ(nodeInLevelToGridPointRef6[i],
                  nodeInLevelToGridPointPtr[i]);
    }

    // Sweep direction 7
    Graph3D<SweepNumber3D::Sweep7> graph7;
    EXPECT_NO_THROW(graph7.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph7.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {   
        EXPECT_EQ(nodeInLevelToGridPointRef7[i],
                  nodeInLevelToGridPointPtr[i]);
    }

    // Sweep direction 8
    Graph3D<SweepNumber3D::Sweep8> graph8;
    EXPECT_NO_THROW(graph8.initialize(nx, ny, nz));
    nodeInLevelToGridPointPtr = graph8.getNodeInLevelToIndexPointer();
    for (int i = 0; i < graph1.getNumberOfGridPoints(); ++i)
    {
        EXPECT_EQ(nodeInLevelToGridPointRef8[i],
                  nodeInLevelToGridPointPtr[i]);
    }
}
}
