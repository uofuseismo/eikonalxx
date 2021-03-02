#ifndef PRIVATE_GRID_HPP
#define PRIVATE_GRID_HPP
namespace
{

/// @brief Converts a grid index (ix,iz) to an index.
/// @param[in] nx   The number of grid points in x.
/// @param[in] ix   The ix'th grid point.
/// @param[in] iz   The iz'th grid point.
/// @result The index in a row-major [nz x nx] matrix corresponding to the
///         given (ix, iz)'th grid point.
#pragma omp declare(simd) uniform(nx)
[[nodiscard]] [[maybe_unused]]
int gridToIndex(const int nx, const int ix, const int iz)
{
    return iz*nx + ix;
}
/// @brief Converts a grid index (ix,iy,iz) to an index.
/// @param[in] nx   The number of grid points in x.
/// @param[in] ny   The number of grid points in y.
/// @param[in] nz   The number of grid points in z.
/// @param[in] ix   The ix'th grid point.
/// @param[in] iy   The iy'th grid point.
/// @param[in] iz   The iz'th grid point
#pragma omp declare(simd) uniform(nx, ny)
[[nodiscard]]
int gridToIndex(const int nx, const int ny,
                const int ix, const int iy, const int iz)
{
    return iz*(nx*ny) + iy*nx + ix;
}

/// @brief This is the inverse operation of gridToIndex and converts an
///        a grid index to an (ix, iy, iz).
/// @param[in] igrd  The grid index.
/// @param[in] nx    The number of grid points in x.
/// @param[in] ny    The number of grid points in y.
/// @param[out] ix   The ix'th grid point.
/// @param[out] iy   The iy'th grid point.
/// @param[out] iz   The iz'th grid point.
#pragma omp declare(simd) uniform(nx, ny)
void indexToGrid(const int igrd,
                 const int nx, const int ny,
                 int *ix, int *iy, int *iz)
{
    auto nxy = nx*ny;
    *iz = igrd/nxy;
    *iy = (igrd - *iz*nxy)/nx;
    *ix = igrd - *iz*nxy - *iy*nx;
}

}
#endif
