#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <filesystem>
#include <cstdint>
#include "eikonalxx/io/vtkRectilinearGrid2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "pack.hpp"
#include "private/grid.hpp"

using namespace EikonalXX::IO;

class VTKRectilinearGrid2D::VTKRectilinearGrid2DImpl
{
public:
    std::fstream mFile;
    EikonalXX::Geometry2D mGeometry;
    bool mHaveGeometry{false};
    bool mWriteBinary{true};
    bool mWroteLookupTable{false};
    bool mWrotePointData{false};
};

/// C'tor
VTKRectilinearGrid2D::VTKRectilinearGrid2D() :
    pImpl(std::make_unique<VTKRectilinearGrid2DImpl> ())
{
}

/// Destructor
VTKRectilinearGrid2D::~VTKRectilinearGrid2D() = default;

/// Open file
void VTKRectilinearGrid2D::open(const std::string &fileName,
                                const EikonalXX::Geometry2D &geometry,
                                const std::string &title,
                                const bool writeBinary)
{
    // If the file was previously open then close it
    close();
    // The following will throw on the geometry
    pImpl->mGeometry = geometry;
    auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
    int ny = 1;
    auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ(); 
    auto dx = static_cast<float> (pImpl->mGeometry.getGridSpacingInX());
    auto dz = static_cast<float> (pImpl->mGeometry.getGridSpacingInZ());
    auto x0 = static_cast<float> (pImpl->mGeometry.getOriginInX()); 
    auto y0 = 0.f;
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
        pImpl->mFile << ::fillBlanksWithUnderscores(title) << std::endl;
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
        auto cx = ::pack(xCoords);
        pImpl->mFile.write(cx.data(), cx.size());
        //auto xPtr = reinterpret_cast<const char *> (xCoords.data());
        //pImpl->mFile.write(xPtr, xCoords.size()*sizeof(float)); 
        pImpl->mFile << std::endl;

        pImpl->mFile << "Y_COORDINATES " << 1 << " float" << std::endl;
        std::vector<float> yCoords(1, 0);
        auto cy = ::pack(yCoords);
        pImpl->mFile.write(cy.data(), cy.size());
        //auto yPtr = reinterpret_cast<const char *> (&y0);
        //pImpl->mFile.write(yPtr, 1*sizeof(float));
        pImpl->mFile << std::endl;

        pImpl->mFile << "Z_COORDINATES " << nz << " float" << std::endl;
        std::vector<float> zCoords(nz);
        for (int i = 0; i < nz; ++i){zCoords[i] =-(z0 + dz*i);}
        auto cz = ::pack(zCoords);
        pImpl->mFile.write(cz.data(), cz.size());
        //auto zPtr = reinterpret_cast<const char *> (zCoords.data());
        //pImpl->mFile.write(zPtr, zCoords.size()*sizeof(float));
        pImpl->mFile << std::endl;
    }
    else
    {
        pImpl->mFile << "X_COORDINATES " << nx << " float" << std::endl;
        for (int ix = 0; ix < nx; ++ix)
        {
            pImpl->mFile << x0 + dx*ix << std::endl;
        }

        pImpl->mFile << "Y_COORDINATES " << 1 << " float" << std::endl;
        pImpl->mFile << y0 << std::endl;

        pImpl->mFile << "Z_COORDINATES " << nz << " float" << std::endl;
        for (int iz = 0; iz < nz; ++iz)
        {
            pImpl->mFile << z0 - dz*iz << std::endl;
        }
    }
}

/// File is open?
bool VTKRectilinearGrid2D::isOpen() const noexcept
{
    return pImpl->mFile.is_open();
}

/// Close the file
void VTKRectilinearGrid2D::close() noexcept
{
    if (isOpen()){pImpl->mFile.close();}
    pImpl->mGeometry.clear();
    pImpl->mHaveGeometry = false;
    pImpl->mWriteBinary = true;
    pImpl->mWrotePointData = false;
    pImpl->mWroteLookupTable = false;
}

/// Write float nodal dataset
template<>
void VTKRectilinearGrid2D::writeNodalDataset(
    const std::string &name,
    const float *data,
    const Ordering2D ordering) const
{
    if (!isOpen())
    {
        std::cerr << "File not open for writing" << std::endl;
    }   
    if (data == nullptr){throw std::invalid_argument("data is NULL");}
    auto nGrid = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPoints());
    if (!pImpl->mWrotePointData)
    {   
        pImpl->mFile << "POINT_DATA " << nGrid << std::endl;
        pImpl->mWrotePointData = true;
    }
    constexpr int nComp = 1;
    if (name.empty())
    {
        pImpl->mFile << "SCALARS dataset float " << nComp << std::endl;
    }
    else
    {
        pImpl->mFile << "SCALARS " << ::fillBlanksWithUnderscores(name)
                     << " float " << nComp << std::endl;
    }
    if (!pImpl->mWroteLookupTable)
    {
        pImpl->mFile << "LOOKUP_TABLE default" << std::endl;
        pImpl->mWroteLookupTable = true;
    }
    if (pImpl->mWriteBinary)
    {
        std::vector<char> cData(nGrid*4);
        if (ordering == EikonalXX::Ordering2D::Natural)
        {
            ::pack(nGrid, data, cData.data());
        }
        else
        {
            auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
            auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
            for (int ix = 0; ix < nx; ++ix)
            {
                for (int iz = 0; iz < nz; ++iz)
                {
                    auto isrc = ::gridToIndex(nx, ix, iz);
                    auto idst = 4*::gridToIndex(nz, iz, ix);
                    ::pack(data[isrc], &cData[idst]);
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

/// Write float nodal vector dataset
template<typename T>
void VTKRectilinearGrid2D::writeNodalVectorDataset(
    const std::string &name,
    const T *data,
    const Ordering2D ordering) const
{
    if (!isOpen())
    {   
        std::cerr << "File not open for writing" << std::endl;
    }   
    if (data == nullptr){throw std::invalid_argument("data is NULL");}
    auto nGrid = static_cast<size_t> (pImpl->mGeometry.getNumberOfGridPoints());
    if (!pImpl->mWrotePointData)
    {
        pImpl->mFile << "POINT_DATA " << nGrid << std::endl;
        pImpl->mWrotePointData = true;
    }
    if (name.empty())
    {   
        pImpl->mFile << "VECTORS vector_dataset float " << std::endl;
    }
    else
    {
        pImpl->mFile << "VECTORS " << ::fillBlanksWithUnderscores(name)
                     << " float " << std::endl;
    }
    if (!pImpl->mWroteLookupTable)
    {
        pImpl->mFile << "LOOKUP_TABLE default" << std::endl;
        pImpl->mWroteLookupTable = true;
    }
    if (pImpl->mWriteBinary)
    {
        std::vector<char> cData(nGrid*4*3);
        if (ordering == EikonalXX::Ordering2D::Natural)
        {
            auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
            auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
            for (int iz = 0; iz < nz; ++iz)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    auto iSrc = ::gridToIndex(nx, ix, iz);
                    ::pack(data[2*iSrc],     &cData[4*3*iSrc]);
                    ::pack(0.0f,             &cData[4*3*iSrc + 4]); 
                    ::pack(data[2*iSrc + 1], &cData[4*3*iSrc + 8]);
                }   
            }
        }
        else
        {
            auto nx = pImpl->mGeometry.getNumberOfGridPointsInX();
            auto nz = pImpl->mGeometry.getNumberOfGridPointsInZ();
            for (int ix = 0; ix < nx; ++ix)
            {
                for (int iz = 0; iz < nz; ++iz)
                {
                    auto iSrc = ::gridToIndex(nx, ix, iz);
                    auto iDst = ::gridToIndex(nz, iz, ix);
                    ::pack(data[2*iSrc],     &cData[4*3*iDst]);
                    ::pack(0.0f,             &cData[4*3*iDst + 4]);
                    ::pack(data[2*iSrc + 1], &cData[4*3*iDst + 8]);
                }
            }
        }
        pImpl->mFile.write(cData.data(), cData.size());
    }
    else
    {
        for (int iGrd = 0; iGrd < static_cast<int> (nGrid); ++iGrd)
        {
            pImpl->mFile << data[2*iGrd] << " 0 "
                         << data[2*iGrd + 1] << std::endl;
        }
    }
}

/// Write float cellular dataset
template<>
void VTKRectilinearGrid2D::writeCellularDataset(
    const std::string &name,
    const float *data,
    const Ordering2D ordering) const
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
        pImpl->mFile << "SCALARS " << ::fillBlanksWithUnderscores(name)
                     << " float " << nComp << std::endl;
    }
    if (!pImpl->mWroteLookupTable)
    {
        pImpl->mFile << "LOOKUP_TABLE default" << std::endl;
        pImpl->mWroteLookupTable = true;
    }
    if (pImpl->mWriteBinary)
    { 
        std::vector<char> cData(nCell*4);
        if (ordering == EikonalXX::Ordering2D::Natural)
        {
            ::pack(nCell, data, cData.data());
        }
        else
        {
            auto nCellX = pImpl->mGeometry.getNumberOfCellsInX();
            auto nCellZ = pImpl->mGeometry.getNumberOfCellsInZ();
            for (int ix = 0; ix < nCellX; ++ix)
            {
                for (int iz = 0; iz < nCellZ; ++iz)
                {
                    auto isrc = ::gridToIndex(nCellX, ix, iz);
                    auto idst = 4*::gridToIndex(nCellZ, iz, ix);
                    ::pack(data[isrc], &cData[idst]);
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
void VTKRectilinearGrid2D::writeNodalDataset(
    const std::string &name,
    const double *__restrict__ data,
    const Ordering2D ordering) const
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
void VTKRectilinearGrid2D::writeCellularDataset(
    const std::string &name,
    const double *__restrict__ data,
    const Ordering2D ordering) const
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

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template void EikonalXX::IO::VTKRectilinearGrid2D::writeNodalVectorDataset(
    const std::string &, const float *, EikonalXX::Ordering2D) const;
template void EikonalXX::IO::VTKRectilinearGrid2D::writeNodalVectorDataset(
    const std::string &, const double *, EikonalXX::Ordering2D) const;

