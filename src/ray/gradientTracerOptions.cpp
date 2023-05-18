#include <vector>
#include <algorithm>
#include <limits>
#include "eikonalxx/ray/gradientTracerOptions.hpp"

using namespace EikonalXX::Ray;

class GradientTracerOptions::GradientTracerOptionsImpl
{
public:
    std::vector<std::pair<int, double>> mRadiusScaleFactor
    {
        std::pair{1, 0.05},
        std::pair{2, 0.10},
        std::pair{3, 0.20},
        std::pair{5, 0.40},
        std::pair{std::numeric_limits<int>::max(), 0.45}
    };
};

/// Constructor
GradientTracerOptions::GradientTracerOptions() :
    pImpl(std::make_unique<GradientTracerOptionsImpl> ())
{
}

/// Copy constructor
GradientTracerOptions::GradientTracerOptions(
    const GradientTracerOptions &options)
{
    *this = options;
}

/// Move constructor
GradientTracerOptions::GradientTracerOptions(
    GradientTracerOptions &&options) noexcept
{   
    *this = std::move(options);
}

/// Copy assignment
GradientTracerOptions&
GradientTracerOptions::operator=(const GradientTracerOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<GradientTracerOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
GradientTracerOptions&
GradientTracerOptions::operator=(GradientTracerOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Reset class
void GradientTracerOptions::clear() noexcept
{
    pImpl = std::make_unique<GradientTracerOptionsImpl> ();
}

/// Destructor
GradientTracerOptions::~GradientTracerOptions() = default;

void GradientTracerOptions::setRadiusScaleFactor(
    const std::vector<std::pair<int, double>> &radiusScaleFactor)
{
    if (radiusScaleFactor.empty())
    {
        throw std::invalid_argument("Input vector cannot be empty");
    }
    for (const auto &item : radiusScaleFactor)
    {
        if (item.first < 0)
        {
            throw std::invalid_argument("All radii must be non-negative");
        }
        if (item.second <= std::numeric_limits<double>::epsilon())
        {
            throw std::invalid_argument(
                "All scale factors must be non-negative");
        }
        if (item.second >= 1)
        {
            throw std::invalid_argument(
                "All scale factors must be less than 1");
        }
    }
    auto temp = radiusScaleFactor;
    std::sort(temp.begin(), temp.end(),
              [](const std::pair<int, double> &lhs,
                 const std::pair<int, double> &rhs)
              {
                  return lhs.first < rhs.first;
              });
    for (int i = 0; i < static_cast<int> (temp.size() - 1); ++i)
    {
        if (temp[i].first == temp[i+1].first)
        {
            throw std::invalid_argument("Duplicate radii");
        }
    }
    if (temp.back().first < std::numeric_limits<int>::max())
    {
        temp.push_back( std::pair{std::numeric_limits<int>::max(),
                                  temp.back().second});
    }
    pImpl->mRadiusScaleFactor = temp; 
}

std::vector<std::pair<int, double>>
GradientTracerOptions::getRadiusScaleFactor() const noexcept
{
    return pImpl->mRadiusScaleFactor;
}
