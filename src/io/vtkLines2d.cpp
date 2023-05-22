#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include "eikonalxx/io/vtkLines2d.hpp"
#include "eikonalxx/geometry2d.hpp"
#include "pack.hpp"

using namespace EikonalXX::IO;

class VTKLines2D::VTKLines2DImpl
{
public:
    std::fstream mFile;
    EikonalXX::Geometry2D mGeometry;
    bool mWriteBinary{false};
};

/// Constructor
VTKLines2D::VTKLines2D() :
    pImpl(std::make_unique<VTKLines2DImpl> ())
{
}

/// Destructor
VTKLines2D::~VTKLines2D() = default;

/// Open file
void VTKLines2D::open(const std::string &fileName,
                        const EikonalXX::Geometry2D &geometry,
                        const std::string &title)
{
    constexpr bool writeBinary{false};
    // If the file was previously open then close it
    close();
    pImpl->mGeometry = geometry;
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
void VTKLines2D::write(
    const std::vector<std::vector<std::pair<double, double>>> &polygons) const
{
    if (!isOpen()){throw std::runtime_error("File not open");}
    int nPoints = 0;
    for (const auto &p : polygons)
    {
        nPoints = nPoints + static_cast<int> (p.size());
    }
    pImpl->mFile << "POINTS " << nPoints << " float" << std::endl;
    auto z0 = pImpl->mGeometry.getOriginInZ();
    for (const auto &p : polygons)
    {
        for (const auto &pi : p)
        {
            auto zi =-z0 - (pi.second - z0); // Reverse
            pImpl->mFile << pi.first << " 0 " << zi << std::endl;
        }
    }
    pImpl->mFile << "LINES " << polygons.size() << " "
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
bool VTKLines2D::isOpen() const noexcept
{
    return pImpl->mFile.is_open();
}

/// Close the file
void VTKLines2D::close() noexcept
{
    if (isOpen()){pImpl->mFile.close();}
}
