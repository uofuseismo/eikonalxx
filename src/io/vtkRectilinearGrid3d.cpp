#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <filesystem>
#include <cstdint>
#include "eikonalxx/io/vtkRectilinearGrid3d.hpp"
#include "eikonalxx/geometry3d.hpp"
#include "private/grid.hpp"
#include "private/pack.hpp"

using namespace EikonalXX::IO;

class VTKRectilinearGrid3D::VTKRectilinearGrid3DImpl
{
public:
    std::fstream mFile;
    EikonalXX::Geometry3D mGeometry;
    bool mHaveGeometry = false;
    bool mWriteBinary = true;
};

/// C'tor
VTKRectilinearGrid3D::VTKRectilinearGrid3D() :
    pImpl(std::make_unique<VTKRectilinearGrid3DImpl> ())
{
}

/// Destructor
VTKRectilinearGrid3D::~VTKRectilinearGrid3D() = default;

/// Open file
void VTKRectilinearGrid3D::open(const std::string &fileName,
                                const EikonalXX::Geometry3D &geometry,
                                const std::string &title,
                                const bool writeBinary)
{
    // If the file was previously open then close it
    close();
    // The following will throw on the geometry
    pImpl->mGeometry = geometry;
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    auto ny = pImpl->mGeometry.getNumberOfGridPointsInY();
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ(); 
    auto dx = static_cast<float> (pImpl->mGeometry.getGridSpacingInX());
    auto dy = static_cast<float> (pImpl->mGeometry.getGridSpacingInY());
    auto dz = static_cast<float> (pImpl->mGeometry.getGridSpacingInZ());
    auto x0 = static_cast<float> (pImpl->mGeometry.getOriginInX()); 
    auto y0 = static_cast<float> (pImpl->mGeometry.getOriginInY());
    auto z0 = static_cast<float> (pImpl->mGeometry.getOriginInZ()); 
    // Warn user file will be overwritten if it exists
    if (std::filesystem::exists(fileName))
    {
        std::cerr << "Overwriting: " << fileName << std::endl;
    }
    else
    {
        // If directory containing file doesn't exist then make it
        auto parentPath = std::filesystem::path(fileName).parent_path();
        if (!parentPath.empty() && !std::filesystem::exists(parentPath))
        {
            std::cerr << "Making path: " << parentPath << std::endl;
            std::filesystem::create_directory(parentPath);
        }
    }
    // Open file in binary or ASCII
    pImpl->mWriteBinary = writeBinary;
    if (pImpl->mWriteBinary)
    {
        pImpl->mFile = std::fstream(fileName, std::ios::out | std::ios::trunc |
                                              std::ios::binary);
    }
    else
    {
        pImpl->mFile = std::fstream(fileName, std::ios::out | std::ios::trunc);
    }
    // (1) Header
    pImpl->mFile << "# vtk DataFile Version 3.0" << std::endl;
    // (2) Title
    if (title.empty())
    {
        pImpl->mFile << "vtk_dataset" << std::endl;
    }
    else
    {
        pImpl->mFile << fillBlanksWithUnderscores(title) << std::endl;
    }
    // (3) Data type
    if (pImpl->mWriteBinary)
    {
        pImpl->mFile << "BINARY" << std::endl;
    }
    else
    {
        pImpl->mFile << "ASCII" << std::endl;
    } 
    // (4) Dataset (geometry)
    pImpl->mFile << "DATASET RECTILINEAR_GRID" << std::endl;
    pImpl->mFile << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
    if (pImpl->mWriteBinary)
    {
        pImpl->mFile << "X_COORDINATES " << nx << " float" << std::endl;
        std::vector<float> xCoords(nx); 
        for (int i = 0; i < nx; ++i){xCoords[i] = x0 + dx*i;}
        auto cx = pack(xCoords);
        pImpl->mFile.write(cx.data(), cx.size());
        pImpl->mFile << std::endl;

        pImpl->mFile << "Y_COORDINATES " << ny << " float" << std::endl;
        std::vector<float> yCoords(ny);
        for (int i = 0; i < ny; ++i){yCoords[i] = y0 + dy*i;}
        auto cy = pack(yCoords);
        pImpl->mFile.write(cy.data(), cy.size());
        pImpl->mFile << std::endl;

        pImpl->mFile << "Z_COORDINATES " << nz << " float" << std::endl;
        std::vector<float> zCoords(nz);
        for (int i = 0; i < nz; ++i){zCoords[i] = z0 - dz*i;}
        auto cz = pack(zCoords);
        pImpl->mFile.write(cz.data(), cz.size());
        pImpl->mFile << std::endl;
    }
    else
    {
        pImpl->mFile << "X_COORDINATES " << nx << " float" << std::endl;
        for (int ix = 0; ix < nx; ++ix)
        {
            pImpl->mFile << x0 + dx*ix << std::endl;
        }

        pImpl->mFile << "Y_COORDINATES " << ny << " float" << std::endl;
        for (int iy = 0; iy < ny; ++iy)
        {
            pImpl->mFile << y0 + dy*iy << std::endl;
        }

        pImpl->mFile << "Z_COORDINATES " << nz << " float" << std::endl;
        for (int iz = 0; iz < nz; ++iz)
        {
            pImpl->mFile << z0 - dz*iz << std::endl;
        }
    }
}

/// File is open?
bool VTKRectilinearGrid3D::isOpen() const noexcept
{
    return pImpl->mFile.is_open();
}

/// Close the file
void VTKRectilinearGrid3D::close() noexcept
{
    if (isOpen()){pImpl->mFile.close();}
    pImpl->mGeometry.clear();
    pImpl->mHaveGeometry = false;
    pImpl->mWriteBinary = true;
}

/// Write float nodal dataset
template<>
void VTKRectilinearGrid3D::writeNodalDataset(
    const std::string &name,
    const float *data,
    const Ordering3D ordering) const
{
    if (!isOpen())
    {
        std::cerr << "File not open for writing" << std::endl;
    }   
    if (data == nullptr){throw std::invalid_argument("data is NULL");}
    auto nGrid = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPoints());
    pImpl->mFile << "POINT_DATA " << nGrid << std::endl;
    constexpr int nComp = 1;
    if (name.empty())
    {
        pImpl->mFile << "SCALARS dataset float " << nComp << std::endl;
    }
    else
    {
        pImpl->mFile << "SCALARS " << fillBlanksWithUnderscores(name)
                     << " float " << nComp << std::endl;
    }
    pImpl->mFile << "LOOKUP_TABLE default" << std::endl;
    if (pImpl->mWriteBinary)
    {
        std::vector<char> cData(nGrid*4);
        if (ordering == EikonalXX::Ordering3D::Natural)
        {
            pack(nGrid, data, cData.data());
        }
        else
        {
            auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
            auto ny = pImpl->mGeometry.getNumberOfGridPointsInY();
            auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
            for (int ix = 0; ix < nx; ++ix)
            {
                for (int iy = 0; iy < ny; ++iy)
                {
                    for (int iz = 0; iz < nz; ++iz)
                    {
                        auto isrc = gridToIndex(nx, ny, ix, iy, iz);
                        auto idst = 4*gridToIndex(nz, ny, iz, iy, ix);
                        pack(data[isrc], &cData[idst]);
                    }
                }
            }
        }
        pImpl->mFile.write(cData.data(), cData.size());
    }
    else
    {
        for (size_t i = 0; i < nGrid; ++i)
        {
            pImpl->mFile << data[i] << std::endl;
        }
    }
}

/// Write float nodal dataset
template<>
void VTKRectilinearGrid3D::writeCellularDataset(
    const std::string &name,
    const float *data,
    const Ordering3D ordering) const
{
    if (!isOpen())
    {
        std::cerr << "File not open for writing" << std::endl;
    }
    if (data == nullptr){throw std::invalid_argument("data is NULL");}
    auto nCell = static_cast<size_t> (pImpl->mGeometry.getNumberOfCells());
    pImpl->mFile << "CELL_DATA " << nCell << std::endl;
    constexpr int nComp = 1;
    if (name.empty())
    {
        pImpl->mFile << "SCALARS dataset float " << nComp << std::endl;
    }
    else
    {
        pImpl->mFile << "SCALARS " << fillBlanksWithUnderscores(name)
                     << " float " << nComp << std::endl;
    }
    pImpl->mFile << "LOOKUP_TABLE default" << std::endl;
    if (pImpl->mWriteBinary)
    { 
        std::vector<char> cData(nCell*4);
        if (ordering == EikonalXX::Ordering3D::Natural)
        {
            pack(nCell, data, cData.data());
        }
        else
        {
            auto nCellX = pImpl->mGeometry.getNumberOfCellsInX();
            auto nCellY = pImpl->mGeometry.getNumberOfCellsInY();
            auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
            for (int ix = 0; ix < nCellX; ++ix)
            {
                for (int iy = 0; iy < nCellY; ++iy)
                {
                    for (int iz = 0; iz < nCellZ; ++iz)
                    {
                        auto isrc = gridToIndex(nCellX, nCellY, ix, iy, iz);
                        auto idst = 4*gridToIndex(nCellZ, nCellY, iz, iy, ix);
                        pack(data[isrc], &cData[idst]);
                    }
                }
            }
        }
        pImpl->mFile.write(cData.data(), cData.size());
    }
    else
    {
        for (size_t i = 0; i < nCell; ++i)
        {
            pImpl->mFile << data[i] << std::endl;
        }
    }
}

/// Double nodal dataset
template<>
void VTKRectilinearGrid3D::writeNodalDataset(
    const std::string &name,
    const double *__restrict__ data,
    const Ordering3D ordering) const
{
    if (!isOpen())
    {   
        std::cerr << "File not open for writing" << std::endl;
    }   
    auto nGrid = pImpl->mGeometry.getNumberOfGridPoints();
    std::vector<float> data4(nGrid);
    float *__restrict__ data4Ptr = data4.data();
    std::copy(data, data + nGrid, data4Ptr);
    writeNodalDataset(name, data4Ptr, ordering);
}

/// Double cell data
template<>
void VTKRectilinearGrid3D::writeCellularDataset(
    const std::string &name,
    const double *__restrict__ data,
    const Ordering3D ordering) const
{
    if (!isOpen())
    {
        std::cerr << "File not open for writing" << std::endl;
    }
    auto nCell = pImpl->mGeometry.getNumberOfCells();
    std::vector<float> data4(nCell);
    float *__restrict__ data4Ptr = data4.data();
    std::copy(data, data + nCell, data4Ptr);
    writeCellularDataset(name, data4Ptr, ordering);
}
