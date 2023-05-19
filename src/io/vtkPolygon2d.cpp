#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include "eikonalxx/io/vtkPolygon2d.hpp"
#include "pack.hpp"

using namespace EikonalXX::IO;

class VTKPolygon2D::VTKPolygon2DImpl
{
public:
    std::fstream mFile;
    bool mWriteBinary{false};
};

/// Constructor
VTKPolygon2D::VTKPolygon2D() :
    pImpl(std::make_unique<VTKPolygon2DImpl> ())
{
}

/// Destructor
VTKPolygon2D::~VTKPolygon2D() = default;

/// Open file
void VTKPolygon2D::open(const std::string &fileName,
                        const std::string &title)
{
    constexpr bool writeBinary{false};
    // If the file was previously open then close it
    close();
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
    pImpl->mFile << "DATASET POLYDATA " << std::endl;
}

/// Write the polygons
void VTKPolygon2D::write(
    const std::vector<std::vector<std::pair<double, double>>> &polygons) const
{
    if (!isOpen()){throw std::runtime_error("File not open");}
    int nPoints = 0;
    for (const auto &p : polygons)
    {
        nPoints = nPoints + static_cast<int> (p.size());
    }
    pImpl->mFile << "POINTS " << nPoints << " float" << std::endl;
    for (const auto &p : polygons)
    {
        for (const auto &pi : p)
        {
            pImpl->mFile << pi.first << " 0 " << pi.second << std::endl;
        }
    }
    pImpl->mFile << "POLYGONS " << polygons.size() << " "
                                << nPoints + polygons.size() << std::endl; 
    int iPoint = 0;
    for (const auto &p : polygons)
    {
        pImpl->mFile << p.size() << " ";//"4 ";
        for (int j = 0; j < static_cast<int> (p.size()); ++j)
        { 
            pImpl->mFile << iPoint << " ";
            iPoint = iPoint + 1;
        }
        pImpl->mFile << std::endl;
    }
}

/// File is open?
bool VTKPolygon2D::isOpen() const noexcept
{
    return pImpl->mFile.is_open();
}

/// Close the file
void VTKPolygon2D::close() noexcept
{
    if (isOpen()){pImpl->mFile.close();}
}
