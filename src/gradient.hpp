#ifndef PRIVATE_GRADIENT_2D_HPP
#define PRIVATE_GRADIENT_2D_HPP
#include <cmath>
#include <CL/sycl.hpp>
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
#include "eikonalxx/model2d.hpp"
#include "private/grid.hpp"

namespace
{

enum class DerivativeType
{
    ForwardDifference = 0,
    CentralDifference = 1
};

template<typename T>
void finiteDifference(const EikonalXX::Geometry2D &geometry,
                      const EikonalXX::Source2D &source,
                      const EikonalXX::Model2D<T> &velocityModel,
                      const T *travelTimeField,
                      std::vector<T> *gradientInX,
                      std::vector<T> *gradientInZ,
                      const DerivativeType derivativeType = DerivativeType::CentralDifference)
{
    // Basic size check
    auto nGridX = static_cast<size_t> (geometry.getNumberOfGridPointsInX());
    auto nGridZ = static_cast<size_t> (geometry.getNumberOfGridPointsInZ());
    if (derivativeType == DerivativeType::CentralDifference)
    {
        if (nGridX < 3)
        {
            throw std::invalid_argument(
                "At least 3 grid points in x for central difference");
        }
        if (nGridZ < 3)
        {
            throw std::invalid_argument(
                "At least 3 grid points in z for central difference");
        }
    }
    else
    {
        if (nGridX < 2)
        {
            throw std::invalid_argument(
                "At least 3 grid points in x for forward/backward difference");
        }
        if (nGridZ < 2)
        {
            throw std::invalid_argument(
                "At least 2 grid points in z for forward/backward difference");
        }
    }
    // Queue
    sycl::queue q{sycl::cpu_selector_v,
                  sycl::property::queue::in_order()};
    auto workGroupSize = q.get_device().get_info<sycl::info::device::max_work_group_size> (); 
    workGroupSize = static_cast<size_t> (std::sqrt(workGroupSize));
    // Space allocation
    auto nGrid = nGridX*nGridZ;
    if (gradientInX->size() != nGrid)
    {
        gradientInX->resize(nGrid, 0);
    }
    if (gradientInZ->size() != nGrid)
    {
        gradientInZ->resize(nGrid, 0);
    }
    auto dx = geometry.getGridSpacingInX();
    auto dz = geometry.getGridSpacingInZ();
    auto dxi = static_cast<T> (static_cast<double> (1./dx));
    auto dzi = static_cast<T> (static_cast<double> (1./dz));
    auto twodxi = static_cast<T> (1./(2.*dx));
    auto twodzi = static_cast<T> (1./(2.*dz));
    // Determine ranges
    sycl::range global{nGridX, nGridZ};
    auto nTileX = sycl::min(workGroupSize, nGridX);
    auto nTileZ = sycl::min(workGroupSize, nGridZ);
    sycl::range local{nTileX, nTileZ};
    // Create the buffers
    {
    sycl::buffer<T> gradientBufferX(*gradientInX);
    sycl::buffer<T> gradientBufferZ(*gradientInZ);
    sycl::buffer<T> travelTimeFieldBuffer(travelTimeField, nGrid);
    // Tabulate the appropriate gradient which involves taking the appropriate
    // derivative.  Note, first focus on the volume then clean up the edges
    // and source after the fact.
    auto eGradient = q.submit([&](sycl::handler &h) 
    {
        sycl::accessor gradientInXAccessor(gradientBufferX, h,
                                           sycl::write_only, sycl::no_init);
        sycl::accessor gradientInZAccessor(gradientBufferZ, h,
                                           sycl::write_only, sycl::no_init);
        sycl::accessor travelTimeAccessor(travelTimeFieldBuffer, h,
                                          sycl::read_only);
        if (derivativeType == DerivativeType::CentralDifference)
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<2> it)
            {
                size_t ix = it.get_global_id(0);
                size_t iz = it.get_global_id(1);
                auto index = ::gridToIndex(nGridX, ix, iz);

                T uxm1 = 0;
                if (ix > 0){uxm1 = travelTimeAccessor[index - 1];}
                T uxp1 = 0;
                if (ix < nGridX - 1){uxp1 = travelTimeAccessor[index + 1];}
                gradientInXAccessor[index] = (uxp1 - uxm1)*twodxi;

                T uzm1 = 0;
                if (iz > 0){uzm1 = travelTimeAccessor[index - nGridX];}
                T uzp1 = 0;
                if (iz < nGridZ - 1){uzp1 = travelTimeAccessor[index + nGridX];}
                gradientInZAccessor[index] = (uzp1 - uzm1)*twodxi;
            });
        }
        else if (derivativeType == DerivativeType::ForwardDifference)
        {
            h.parallel_for(sycl::nd_range{global, local},
                           [=](sycl::nd_item<2> it) 
            {
                size_t ix = it.get_global_id(0);
                size_t iz = it.get_global_id(1);
                auto index = ::gridToIndex(nGridX, ix, iz);

                T u = travelTimeAccessor[index];
                T uxp1 = 0;
                if (ix < nGridX - 1){uxp1 = travelTimeAccessor[index + 1];}
                gradientInXAccessor[index] = (uxp1 - u)*dxi;

                T uzp1 = 0;
                if (iz < nGridZ - 1){uzp1 = travelTimeAccessor[index + nGridX];}
                gradientInZAccessor[index] = (uzp1 - u)*dzi;
            });
        }
    });
    } // Buffers go out of scope which transfers data back to host
    if (derivativeType == DerivativeType::CentralDifference)
    {
        // Clean it up - backward and forward difference.  The Lagrange
        // polynomial is:
        // P_2(x) = u_0 (x - h)*(x - 2h)/(-2h^2)
        //        + u_1 x*(x - 2h)/(-h^2)
        //        + u_2 x*(x - h)/(2h^2)
        // Differentiating w.r.t. x yields:
        //   P_2(x)' = ( u_0*(2x - 3h) + 4*u_1*(h - x) - u_2 (h - 2x) )/(2h^2)
        // Evaluate at 0 to obtain
        //   dP_2/dx |x=0  = -3 u_0/(2h) + 4 u_1/(2h) - u_2/(2h) 
        // Evaluate at 2h to obtain
        //   dP_2/dx |x=2h = -u_0/(2h) - 4 u_1/(2h) + 3 u_2/(2h)
        // This is the formula given in Mathews - Numerical Methods for Math,
        // Science, and Engineering
        T *gradientInXPtr = gradientInX->data(); 
        for (size_t iz = 0; iz < nGridZ; ++iz)
        {
            constexpr size_t ix = 0;
            auto index = ::gridToIndex(nGridX, ix, iz);
            T u = travelTimeField[index]; 
            T uxp1 = travelTimeField[index + 1];
            T uxp2 = travelTimeField[index + 2];
            gradientInXPtr[index] = (-3*u + 4*uxp1 - uxp2)*twodxi;

            index = ::gridToIndex(nGridX, nGridX - 1, iz);
            T uxm2 = travelTimeField[index - 2];
            T uxm1 = travelTimeField[index - 1];
            u = travelTimeField[index];
            gradientInXPtr[index] = ( uxm2 - 4*uxm1 + 3*u)*twodxi;
        }
        T *gradientInZPtr = gradientInZ->data();
        for (size_t ix = 0; ix < nGridX; ++ix)
        {
            constexpr size_t iz = 0;
            auto index = ::gridToIndex(nGridX, ix, iz);
            T u = travelTimeField[index];
            T uzp1 = travelTimeField[index +   nGridX];
            T uzp2 = travelTimeField[index + 2*nGridX];
            gradientInZPtr[index] = (-3*u + 4*uzp1 - uzp2)*twodzi;

            index = ::gridToIndex(nGridX, ix, nGridZ - 1);
            T uzm2 = travelTimeField[index - 2*nGridX];
            T uzm1 = travelTimeField[index -   nGridX];
            u = travelTimeField[index];
            gradientInZPtr[index] = ( uzm2 - 4*uzm1 + 3*u)*twodzi;
        }
    }
    else if (derivativeType == DerivativeType::ForwardDifference)
    {
        T *gradientInXPtr = gradientInX->data();
        for (size_t iz = 0; iz < nGridZ; ++iz)
        {
            constexpr size_t ix = 0;
            auto index = ::gridToIndex(nGridX, ix, iz);
            T u = travelTimeField[index];
            T uxp1 = travelTimeField[index + 1];
            gradientInXPtr[index] = (uxp1 - u)*dxi;

            index = ::gridToIndex(nGridX, nGridX - 1, iz);
            T uxm1 = travelTimeField[index - 1];
            u = travelTimeField[index];
            gradientInXPtr[index] = (u - uxm1)*dxi;
        }
        T *gradientInZPtr = gradientInZ->data();
        for (size_t ix = 0; ix < nGridX; ++ix)
        {
            constexpr size_t iz = 0;
            auto index = ::gridToIndex(nGridX, ix, iz);
            T u = travelTimeField[index];
            T uzp1 = travelTimeField[index + nGridX];
            gradientInZPtr[index] = (uzp1 - u)*dzi;

            index = ::gridToIndex(nGridX, ix, nGridZ - 1);
            T uzm1 = travelTimeField[index - nGridX];
            u = travelTimeField[index];
            gradientInZPtr[index] = (u - uzm1)*dzi;
        }
    }
    // Fix the derivatives around the source with a `forward' difference.
    // The idea is to make this upwinding since there's a discontinuity in the
    // gradient at the source.   
    double xs = source.getOffsetInX(); 
    double zs = source.getOffsetInZ();
    auto sourceCellInX = static_cast<size_t> (source.getCellInX());
    auto sourceCellInZ = static_cast<size_t> (source.getCellInZ());
    auto sourceCell    = static_cast<size_t> (source.getCell());
    // Slowness at cell
    const auto *slowPtr = velocityModel.getSlownessPointer();
    auto slow0 = slowPtr[sourceCell];
    constexpr double tol{1.e-7};
    // Geometric terms
    auto dxs = static_cast<double> (sourceCellInX)*dx - xs;
    auto dzs = static_cast<double> (sourceCellInZ)*dz - zs;
    double dx0 = sourceCellInX*dx - xs; 
    double dx1 = (sourceCellInX + 1)*dx - xs; 
    double dz0 = sourceCellInZ*dz - zs; 
    double dz1 = (sourceCellInZ + 1)*dz - zs;
    //std::cout << sourceCellInX << " " << sourceCellInZ << std::endl;
    // The main thing to fix is the case when the source is in a cell.
    auto i0 = ::gridToIndex(nGridX, sourceCellInX    , sourceCellInZ);
    auto i1 = ::gridToIndex(nGridX, sourceCellInX + 1, sourceCellInZ);
    auto i2 = ::gridToIndex(nGridX, sourceCellInX + 1, sourceCellInZ + 1);
    auto i3 = ::gridToIndex(nGridX, sourceCellInX,     sourceCellInZ + 1);
//std::cout << dxs << " " << dzs << " " << xs << " " << zs << std::endl;
    if (std::abs(dxs) > tol || std::abs(dzs) > tol)
    {
//std::cout << "yar" << std::endl;
        // Analytic gradients assuming constant velocity in cell
        auto gradX00 = (dx0/sycl::sqrt(dx0*dx0 + dz0*dz0))*slow0;
        auto gradZ00 = (dz0/sycl::sqrt(dx0*dx0 + dz0*dz0))*slow0;
        auto gradX10 = (dx1/sycl::sqrt(dx1*dx1 + dz0*dz0))*slow0;
        auto gradZ10 = (dz0/sycl::sqrt(dx1*dx1 + dz0*dz0))*slow0;
        auto gradX11 = (dx1/sycl::sqrt(dx1*dx1 + dz1*dz1))*slow0;
        auto gradZ11 = (dz1/sycl::sqrt(dx1*dx1 + dz1*dz1))*slow0;
        auto gradX01 = (dx0/sycl::sqrt(dx0*dx0 + dz1*dz1))*slow0;
        auto gradZ01 = (dz1/sycl::sqrt(dx0*dx0 + dz1*dz1))*slow0;
        // Fill gradients (cache friendly)
        if (std::abs(dxs) > tol)
        {
            gradientInX->at(i0) = static_cast<T> (gradX00);
            gradientInX->at(i1) = static_cast<T> (gradX10);
            gradientInX->at(i3) = static_cast<T> (gradX01);
            gradientInX->at(i2) = static_cast<T> (gradX11);
        }

        if (std::abs(dzs) > tol)
        {
            gradientInZ->at(i0) = static_cast<T> (gradZ00);
            gradientInZ->at(i1) = static_cast<T> (gradZ10);
            gradientInZ->at(i3) = static_cast<T> (gradZ01);
            gradientInZ->at(i2) = static_cast<T> (gradZ11);
        }
    }
}

}
#endif
