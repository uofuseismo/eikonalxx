#include <iostream>
#include "eikonalxx/solverOptions.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace EikonalXX;

TEST(SolverOptions, solverOptions)
{
    SolverOptions options;
    const double tol{1.e-3};
    const Verbosity verbosity{Verbosity::Info};
    const int nSweeps{12};
    const int epsilon{6};
    const SolverAlgorithm algorithm{SolverAlgorithm::FastSweepingMethod};

    options.setConvergenceTolerance(tol);
    options.setNumberOfSweeps(nSweeps);
    options.setFactoredEikonalEquationSolverRadius(epsilon);
    options.setAlgorithm(algorithm);
    options.setVerbosity(verbosity);
    // Use copy c'tor and check equality
    SolverOptions optionsCopy(options);
    EXPECT_EQ(options.getNumberOfSweeps(), nSweeps);
    EXPECT_EQ(options.getFactoredEikonalEquationSolverRadius(), epsilon);
    EXPECT_NEAR(optionsCopy.getConvergenceTolerance(), tol, 1.e-10);
    EXPECT_EQ(options.getAlgorithm(), algorithm);
    EXPECT_EQ(options.getVerbosity(), verbosity);
}

}
