#ifndef PRIVATE_GRADIENT_2D_HPP
#define PRIVATE_GRADIENT_2D_HPP
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/source2d.hpp"
namespace
{

template<typename T>
void finiteDifference(const EikonalXX::Geometry2D &geometry,
                      const EikonalXX::Source2D &source,
                      const T *travelTimeField,
                      std::vector<T> *gradientInX,
                      std::vector<T> *gradientInZ)
{
    auto nGridX = geometry.getNumberOfGridPointsInX();
    auto nGridZ = geometry.getNumberOfGridPointsInZ();
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
    // Central difference
     
}

}
#endif
