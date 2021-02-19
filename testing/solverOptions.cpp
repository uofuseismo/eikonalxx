#include <iostream>
#include "eikonalxx/solverOptions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(SolverOptions, solverOptions)
{
    SolverOptions options;
    double tol = 1.e-3;
    Verbosity verbosity = Verbosity::INFO;
    int nSweeps = 12;
    int epsilon = 6;
    auto algorithm = SolverAlgorithm::LEVEL_SET_METHOD;

    options.setConvergenceTolerance(tol);
    options.setNumberOfSweeps(nSweeps);
    options.setSphericalSolverRadius(epsilon);
    options.setAlgorithm(algorithm);
    options.setVerbosity(verbosity);
    // Use copy c'tor and check equality
    SolverOptions optionsCopy(options);
    EXPECT_EQ(options.getNumberOfSweeps(), nSweeps);
    EXPECT_EQ(options.getSphericalSolverRadius(), epsilon);
    EXPECT_NEAR(optionsCopy.getConvergenceTolerance(), tol, 1.e-10);
    EXPECT_EQ(options.getAlgorithm(), algorithm);
    EXPECT_EQ(options.getVerbosity(), verbosity);
}

}
