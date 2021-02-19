#ifndef EIKONALXX_ENUMS_HPP
#define EIKONALXX_ENUMS_HPP
namespace EikonalXX
{
/// @brief Defines the ordering of input 2D velocity models.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
enum class Ordering2D
{
    NATURAL = 0, /*! This is the natural ordering in which x is the fastest
                     changing dimension and z is the slowest changing dimension.
                     Necessarily, the leading dimension is nx. */
    ZX = 0,      /*! This is the natural ordering in which x is the fastest
                     changing dimension and z is the slowest changing dimension.
                     Necessarily, the leading dimension is nx. */
    XZ = 1       /*! This is an alternative mdoel ordering in which z is
                     the fastest changing dimension and x is the slowest
                     changing dimension.  In this case, nz is the leading 
                     dimension. */
};

/// @brief Defines the ordering of input 3D velocity models.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
enum class Ordering3D
{
    NATURAL = 0, /*! This is the natural ordering in which x is the fastest,
                     changing dimension, y is the intermediate dimension,
                     and z is the slowest changing dimension.
                     Necessarily, the first leading dimension is nx and ny
                     is the second leading dimension. */
    ZYX = 0,     /*! This is the natural ordering in which x is the fastest,
                     changing dimension, y is the intermediate dimension,
                     and z is the slowest changing dimension.
                     Necessarily, the first leading dimension is nx and ny
                     is the second leading dimension. */
    XYZ = 1      /*! This is an alternative mdoel ordering in which z is
                     the fastest changing dimension, y is the intermediate
                     dimension, and x is the slowest changing dimension.
                     In this case, nz is the first leading dimension and ny
                     is the second leading dimension. */
};

/// @brief Defines the eiknoal solver algorithm.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
enum class SolverAlgorithm
{
    LEVEL_SET_METHOD,     /*! This solves the eikonal equation using a
                              reordering that allows for parallel updates
                              of travel time nodes during a sweep.  */
    FAST_SWEEPING_METHOD  /*! This is the traditional fast-sweeping method
                              that does not provide for parallel updates
                              of travel time nodes during a sweep. */
};

/// @brief Defines the sweep number for the 2D fast-sweeping method solver.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
enum class SweepNumber2D
{
    SWEEP1 = 0,  /*! This sweep increases in x and z. */ 
    SWEEP2 = 1,  /*! This sweep decreases in x and increases in z. */
    SWEEP3 = 2,  /*! This sweep increases in x and decreases in z. */
    SWEEP4 = 3   /*! This sweep decreases in x and z. */
};

/// @brief Defines the sweep number for the 2D fast-sweeping method solver.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
enum class SweepNumber3D
{
    SWEEP1 = 0,  /*! This sweep increases in x, y, and z. */ 
    SWEEP2 = 1,
    SWEEP3 = 2,
    SWEEP4 = 3,
    SWEEP5 = 4,
    SWEEP6 = 5,
    SWEEP7 = 6,
    SWEEP8 = 7   /*! This sweep decreases in x, y, and z. */
};


/*!
 * @brief Defines the verbosity.
 * @copyright Ben Baker (University of Utah) distributed under the MIT license.
 */
enum class Verbosity
{
    ERROR = 0,   /*!< Only errors are displayed. */
    WARNING = 1, /*!< Errors and warnings. */
    INFO = 2,    /*!< General information. */
    DEBUG = 3    /*!< Fine-grained debugging information. */ 
};

}
#endif
