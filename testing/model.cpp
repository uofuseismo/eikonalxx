#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "eikonalxx/model2d.hpp"
#include "eikonalxx/model3d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include <gtest/gtest.h>

namespace
{

template<class T> std::vector<T> invert(const std::vector<T> &x)
{
    T one = 1;
    std::vector<T> xi(x.size());
    for (size_t i = 0; i < xi.size(); ++i)
    {
        xi[i] = one/x[i];
    }
    return xi;
}

template<class T>
std::vector<T> nodalToCell(const int nx, const int nz, 
                           const std::vector<T> &v)
{
    std::vector<T> s((nx - 1)*(nz - 1), 0); 
    const T four = 4;
    for (int iz=0; iz<nz-1; ++iz)
    {
        for (int ix=0; ix<nx-1; ++ix)
        {
            auto i1 = iz*nx + ix;
            auto i2 = i1 + 1;
            auto i3 = i1 + nx;
            auto i4 = i3 + 1;
            auto idst = iz*(nx - 1) + ix;
            T avgVel = (v[i1] + v[i2] + v[i3] + v[i4])/four;
            s[idst] = 1./avgVel;
        }
    }
    return s;
}

template<class T>
std::vector<T> nodalToCell(const int nx, const int ny, const int nz, 
                           const std::vector<T> &v) 
{
    std::vector<T> s((nx - 1)*(ny - 1)*(nz - 1), 0);
    const T eight = 8;
    for (int iz = 0; iz < nz - 1; ++iz)
    {   
        for (int iy = 0; iy < ny - 1; ++iy)
        {
            for (int ix = 0; ix < nx - 1; ++ix)
            {
                auto i1 = iz*nx*ny + iy*nx + ix; 
                auto i2 = i1 + 1;
                auto i3 = i1 + nx;
                auto i4 = i3 + 1;
                auto i5 = i1 + nx*ny;
                auto i6 = i2 + nx*ny;
                auto i7 = i3 + nx*ny;
                auto i8 = i4 + nx*ny;
                auto idst = iz*(nx - 1)*(ny - 1) + iy*(nx - 1) + ix; 
                T avgVel = (v[i1] + v[i2] + v[i3] + v[i4]
                          + v[i5] + v[i6] + v[i7] + v[i8])/eight;
                s[idst] = 1/avgVel;
            }
        }
    }   
    return s;
}


template<typename T>
T infinityNorm(const std::vector<T> &ref, const std::vector<T> &est)
{
    auto y = ref.data();
    auto yhat = est.data();
    auto n = static_cast<int> (ref.size()); 
    T norm = 0;
    #pragma omp simd reduction(max:norm)
    for (int i=0; i<n; ++i)
    {
        norm = std::max(norm, std::abs(y[i] - yhat[i]));
    }
    return norm;
}
                

TEST(TestModel, model2d)
{
    int nx = 420;//410;
    int nz = 220;
    double dx = 1;
    double dz = 1;
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    auto nGrid = geometry.getNumberOfGridPoints();
    auto nCell = geometry.getNumberOfCells(); 
    // Initialize velocity model
    EikonalXX::Model2D<double> model;
    EXPECT_NO_THROW(model.initialize(geometry)); 

    std::vector<double> velocities(nx*nz, 0);
    std::vector<double> velTrans(nx*nz, 0);
    for (int iz=0; iz<nz; ++iz)
    {
        for (int ix=0; ix<nx; ++ix)
        {
            auto idst = iz*nx + ix;
            velocities[idst] = idst + 1;
            auto jdst = ix*nz + iz;
            velTrans[jdst] = idst + 1;
        }
    }
    std::vector<double> vCell((nx - 1)*(nz - 1), 0);
    std::vector<double> vCellTrans((nx - 1)*(nz - 1), 0);
    for (int iz=0; iz<nz-1; ++iz)
    {
        for (int ix=0; ix<nx-1; ++ix)
        {
            auto idst = iz*(nx - 1) + ix;
            auto jdst = ix*(nz - 1) + iz;
            vCell[idst] = idst + 1;
            vCellTrans[jdst] = idst + 1;
        }
    }
    auto sCellRef = invert(vCell);

    // Set nodal model which requires interpolation
    auto slowRef = nodalToCell(nx, nz, velocities);
    EXPECT_NO_THROW(model.setNodalVelocities(nGrid, velocities.data(),
                                             EikonalXX::Ordering2D::Natural));
    EXPECT_TRUE(model.haveVelocities());
    auto slowness = model.getSlowness();
    EXPECT_EQ(slowness.size(), slowRef.size());
    double errmax = infinityNorm(slowRef, slowness);
    EXPECT_LT(errmax, 1.e-10);
    // Try the transpose operator
    EXPECT_NO_THROW(model.setNodalVelocities(nGrid, velTrans.data(),
                                             EikonalXX::Ordering2D::XZ));
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(slowRef, slowness);
    EXPECT_LT(errmax, 1.e-10);

    // Set the cell-based model which requires no interpolation
    EXPECT_NO_THROW(model.setCellularVelocities(nCell, vCell.data(),
                                              EikonalXX::Ordering2D::Natural)); 
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(sCellRef, slowness);
    EXPECT_LT(errmax, 1.e-10);

    EXPECT_NO_THROW(model.setCellularVelocities(nCell, vCellTrans.data(),
                                                EikonalXX::Ordering2D::XZ));
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(sCellRef, slowness);
    EXPECT_LT(errmax, 1.e-10);
}

TEST(TestModel, model3d)
{
    int nx = 50;
    int ny = 51;
    int nz = 60;
    double dx = 1;
    double dy = 1;
    double dz = 1;
    EikonalXX::Geometry3D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInY(ny);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInY(dy);
    geometry.setGridSpacingInZ(dz);
    auto nGrid = geometry.getNumberOfGridPoints();
    auto nCell = geometry.getNumberOfCells(); 
    // Initialize velocity model
    EikonalXX::Model3D<double> model;
    EXPECT_NO_THROW(model.initialize(geometry)); 

    std::vector<double> velocities(nx*ny*nz, 0); 
    std::vector<double> velTrans(nx*ny*nz, 0); 
    for (int iz = 0; iz < nz; ++iz)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                auto idst = iz*nx*ny + iy*nx + ix;
                velocities.at(idst) = idst + 1;
                auto jdst = ix*ny*nz + iy*nz + iz; 
                velTrans.at(jdst) = idst + 1;
            }
        }
    }   
    std::vector<double> vCell((nx - 1)*(ny - 1)*(nz - 1), 0);
    std::vector<double> vCellTrans((nx - 1)*(ny - 1)*(nz - 1), 0); 
    for (int iz = 0; iz < nz - 1; ++iz)
    {
        for (int iy = 0; iy < ny - 1; ++iy)
        {
            for (int ix = 0; ix < nx - 1; ++ix)
            {
                auto idst = iz*(nx - 1)*(ny - 1) + iy*(nx - 1) + ix; 
                auto jdst = ix*(ny - 1)*(nz - 1) + iy*(nz - 1) + iz;
                vCell.at(idst) = idst + 1;
                vCellTrans.at(jdst) = idst + 1;
            }
        }
    }
    auto sCellRef = invert(vCell);

    // Set nodal model which requires interpolation
    auto slowRef = nodalToCell(nx, ny, nz, velocities);
    EXPECT_NO_THROW(model.setNodalVelocities(nGrid, velocities.data(),
                                             EikonalXX::Ordering3D::Natural));
    EXPECT_TRUE(model.haveVelocities());
    auto slowness = model.getSlowness();
    EXPECT_EQ(slowness.size(), slowRef.size());
    double errmax = infinityNorm(slowRef, slowness);
    EXPECT_LT(errmax, 1.e-10);
    // Try the transpose operator
    EXPECT_NO_THROW(model.setNodalVelocities(nGrid, velTrans.data(),
                                             EikonalXX::Ordering3D::XYZ));
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(slowRef, slowness);
    EXPECT_LT(errmax, 1.e-10);

    // Set the cell-based model which requires no interpolation
    EXPECT_NO_THROW(model.setCellularVelocities(nCell, vCell.data(),
                                              EikonalXX::Ordering3D::Natural)); 
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(sCellRef, slowness);
    EXPECT_LT(errmax, 1.e-10);

    EXPECT_NO_THROW(model.setCellularVelocities(nCell, vCellTrans.data(),
                                                EikonalXX::Ordering3D::XYZ));
    EXPECT_TRUE(model.haveVelocities());
    slowness = model.getSlowness();
    errmax = infinityNorm(sCellRef, slowness);
    EXPECT_LT(errmax, 1.e-10);
}

}
