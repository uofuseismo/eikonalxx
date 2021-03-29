#ifndef EIKONALXX_SOURCE3D_HPP
#define EIKONALXX_SOURCE3D_HPP
#include <ostream>
#include <memory>
namespace EikonalXX
{
class Geometry3D;
/// @class Source3D "source3d.hpp" "eikonalxx/source3d.hpp"
/// @brief Defines a source in a 3D model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Source3D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Source3D();
    /// @brief Copy constructor.
    /// @param[in] source  The source location from which to initialize
    ///                    this class.
    Source3D(const Source3D &source);
    /// @brief Move constructor.
    /// @param[in,out] source  The source location from which to initialize
    ///                        this class.  On exit, source's behavior is
    ///                        undefined.
    Source3D(Source3D &&source) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] source  The source location to copy to this.
    /// @result A deep copy of the input source.
    Source3D& operator=(const Source3D &source);
    /// @brief Move assignment operator.
    /// @param[in] source  The source location whose memory will be moved to
    ///                    this.  On exit, source's behavior is undefined.
    /// @result The memory from source moved to this.
    Source3D& operator=(Source3D &&source) noexcept;
    /// @}

    /// @name Initialization
    /// @{
    /// @brief Sets the model geometry.
    /// @param[in] geometry   The model geometry.  At the very least the grid
    ///                       spacing and number of grid points must be set.
    void setGeometry(const Geometry3D &geometry);
    /// @result The geometry.
    [[nodiscard]] Geometry3D getGeometry() const;
    /// @result True indicates that the geometry was set.
    [[nodiscard]] bool haveGeometry() const noexcept;
    /// @}

    /// @name X Location
    /// @{
    /// @param[in] x   The source location in x in meters.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    /// @throws std::invalid_argument if the x position is not in the model.
    void setLocationInX(double x);
    /// @result The source's location in x in meters.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] double getLocationInX() const;
    /// @result The source's offset from the origin in x in meters.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] double getOffsetInX() const;
    /// @result The x cell containing the source.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] int getCellInX() const;
    /// @result True indicates that the x location was set.
    [[nodiscard]] bool haveLocationInX() const noexcept;
    /// @}

    /// @name Y Location
    /// @{
    /// @param[in] y   The source location in y in meters.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    /// @throws std::invalid_argument if the x position is not in the model.
    void setLocationInY(double y);
    /// @result The source's location in y in meters.
    /// @throws std::runtime_error if \c haveLocationInY() is false.
    [[nodiscard]] double getLocationInY() const;
    /// @result The source's offset from the origin in y in meters.
    /// @throws std::runtime_error if \c haveLocationInY() is false.
    [[nodiscard]] double getOffsetInY() const;
    /// @result The y cell containing the source.
    /// @throws std::runtime_error if \c haveLocationInY() is false.
    [[nodiscard]] int getCellInY() const;
    /// @result True indicates that the y location was set.
    [[nodiscard]] bool haveLocationInY() const noexcept;
    /// @}

    /// @name Z Location
    /// @{
    /// @param[in] z    The source location in z in meters.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    /// @throws std::invalid_argument if the z position is not in the model.
    void setLocationInZ(double z);
    /// @brief Sets the z position to the free surface in the rectilinear
    ///        geometry.  This is useful when computing reciprocity travel time
    ///        fields or modeling surface shots.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    void setZToFreeSurface();
    /// @result The source's location in z in meters.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] double getLocationInZ() const;
    /// @result The source's offset from the origin in z in meters.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] double getOffsetInZ() const;
    /// @result The z cell containing the source.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] int getCellInZ() const;
    /// @result True indicates that the z location was set.
    [[nodiscard]] bool haveLocationInZ() const noexcept;
    /// @}

    /// @name Destructor
    /// @{
    /// @brief Destructor.
    ~Source3D();
    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @}
private:
    class Source3DImpl;
    std::unique_ptr<Source3DImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const Source3D &source);
}
#endif
