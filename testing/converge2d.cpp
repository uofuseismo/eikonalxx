#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <array>
#include <eikonalxx/geometry2d.hpp>
#include <eikonalxx/source2d.hpp>
#include <eikonalxx/model2d.hpp>
#include <eikonalxx/solver2d.hpp>
#include <eikonalxx/solverOptions.hpp>
#include <eikonalxx/analytic/homogeneous2d.hpp>
#include <eikonalxx/analytic/linearGradient2d.hpp>

void homogeneousTest(const double dx);

template<typename T>
[[nodiscard]]
T infinityNorm(const std::vector<T> &x, const std::vector<T> &y)
{
    T infinityNorm = 0;
    if (x.size() != y.size())
    {
        std::cerr << "Warning sizes don't match" << std::endl;
    }
    auto n = std::min(x.size(), y.size());
    const T *__restrict__ xPtr = x.data();
    const T *__restrict__ yPtr = y.data();
    for (size_t i = 0; i < n; ++i)
    {
        infinityNorm = std::max(infinityNorm, std::abs(xPtr[i] - yPtr[i]));
    }
    return infinityNorm;
}

template<typename T>
[[nodiscard]]
double normalizedOneNorm(const std::vector<T> &x, const std::vector<T> &y)
{
    double oneNorm = 0;
    if (x.size() != y.size())
    {   
        std::cerr << "Warning sizes don't match" << std::endl;
    }   
    auto n = std::min(x.size(), y.size());
    const T *__restrict__ xPtr = x.data();
    const T *__restrict__ yPtr = y.data();
    for (size_t i = 0; i< n; ++i)
    {
        oneNorm = oneNorm
                + static_cast<double> (std::abs(xPtr[i] - yPtr[i]));
    }
    return oneNorm/n;
}

template<typename T>
[[nodiscard]]
T normalizedTwoNorm(const std::vector<T> &x, const std::vector<T> &y)
{
    double twoNorm = 0;
    if (x.size() != y.size())
    {
        std::cerr << "Warning sizes don't match" << std::endl;
    }
    auto n = std::min(x.size(), y.size());
    const T *__restrict__ xPtr = x.data();
    const T *__restrict__ yPtr = y.data();
    for (size_t i = 0; i< n; ++i)
    {
        auto residual = static_cast<double> (xPtr[i] - yPtr[i]);
        twoNorm = twoNorm + residual*residual;
    }
    return twoNorm/n;
}


int main()
{
    std::array<double, 3> gridSpacing{100, 10, 1};
    std::cout << std::setprecision(16);
    for (const auto dx : gridSpacing)
    {
        homogeneousTest(dx);
    }
    return EXIT_SUCCESS;
}

// Solve on the 10 km x 10 km domain
void homogeneousTest(const double dx)
{
    const double velocity = 1000; // Uniform
    const double xWidth = 1e3;
    const double zWidth = 1e3;
    const double x0 = 0;
    const double z0 = 0;
    const double dz = dx;
    int nx = std::round(xWidth/dx) + 1;
    int nz = std::round(zWidth/dz) + 1;
    double xs = dx/2; //xWidth/2 + 0.25*dx; // Move off node to make it harder for solver
    double zs = zWidth - dx/2; //zWidth/2 + 0.75*dz; // Move off node to make it harder for solver
    EikonalXX::Geometry2D geometry;
    geometry.setNumberOfGridPointsInX(nx);
    geometry.setNumberOfGridPointsInZ(nz);
    geometry.setGridSpacingInX(dx);
    geometry.setGridSpacingInZ(dz);
    geometry.setOriginInX(x0);
    geometry.setOriginInZ(z0);

    std::vector<double> velocities(nx*nz, velocity);
    EikonalXX::Model2D<double> model;
    model.initialize(geometry);
    model.setNodalVelocities(velocities.size(), velocities.data(),
                             EikonalXX::Ordering2D::Natural);

    EikonalXX::Source2D source;
    source.setGeometry(geometry);
    source.setLocationInX(xs);
    source.setLocationInZ(zs);

    // Create a reference solution
    EikonalXX::Analytic::Homogeneous2D<double> reference;
    reference.initialize(geometry);
    reference.setVelocityModel(velocity);
    reference.setSource(source);
    reference.solve();
    auto exactTimes = reference.getTravelTimeField();

    EikonalXX::SolverOptions options;
    options.setVerbosity(EikonalXX::Verbosity::ERROR);
    options.setAlgorithm(EikonalXX::SolverAlgorithm::FastSweepingMethod);
    options.setConvergenceTolerance(0);

    EikonalXX::Solver2D<double> solver;
    solver.initialize(options, geometry);
    solver.setVelocityModel(model);
    solver.setSource(source);
    solver.solve();

    auto estimateTimes = solver.getTravelTimeField();
    //std::cout << estimateTimes[10] << " " << exactTimes[10] << std::endl;
    std::cout << "Grid Spacing, Normalized One Norm, Normalized Two Norm, Infinity Norm: " 
              << dx << " "
              << normalizedOneNorm(exactTimes, estimateTimes) << " " 
              << normalizedTwoNorm(exactTimes, estimateTimes) << " "  
              << infinityNorm(exactTimes, estimateTimes)
              << std::endl;
}
