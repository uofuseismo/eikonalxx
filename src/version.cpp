#include <string>
#include "eikonalxx/version.hpp"

using namespace EikonalXX;

int Version::getMajor() noexcept
{
    return EIKONALXX_MAJOR;
}

int Version::getMinor() noexcept
{
    return EIKONALXX_MINOR;
}

int Version::getPatch() noexcept
{
    return EIKONALXX_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (EIKONALXX_MAJOR < major){return false;}
    if (EIKONALXX_MAJOR > major){return true;}
    if (EIKONALXX_MINOR < minor){return false;}
    if (EIKONALXX_MINOR > minor){return true;}
    if (EIKONALXX_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(std::to_string(getMajor()) + "."
                      + std::to_string(getMinor()) + "."
                      + std::to_string(getPatch()));
    return version;
}
