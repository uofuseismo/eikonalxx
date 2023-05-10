#ifndef EIKONALXX_SOURCE_2D_HPP
#define EIKONALXX_SOURCE_2D_HPP
#include <ostream>
#include <memory>
namespace EikonalXX
{
class Geometry2D;
/// @class Source2D "source2d.hpp" "eikonalxx/source2d.hpp"
/// @brief Defines a source in a 2D model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Source2D
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Source2D();
    /// @brief Copy constructor.
    /// @param[in] source  The source location from which to initialize
    ///                    this class.
    Source2D(const Source2D &source);
    /// @brief Move constructor.
    /// @param[in,out] source  The source location from which to initialize
    ///                        this class.  On exit, source's behavior is
    ///                        undefined.
    Source2D(Source2D &&source) noexcept;
    /// @}

    /// @name Operators
    /// @{
    /// @brief Copy assignment operator.
    /// @param[in] source  The source location to copy to this.
    /// @result A deep copy of the input source.
    Source2D& operator=(const Source2D &source);
    /// @brief Move assignment operator.
    /// @param[in] source  The source location whose memory will be moved to
    ///                    this.  On exit, source's behavior is undefined.
    /// @result The memory from source moved to this.
    Source2D& operator=(Source2D &&source) noexcept;
    /// @}

    /// @name Initialization
    /// @{
    /// @brief Sets the model geometry.
    /// @param[in] geometry   The model geometry.  At the very least the grid
    ///                       spacing and number of grid points must be set.
    void setGeometry(const Geometry2D &geometry);
    /// @result The geometry.
    [[nodiscard]] Geometry2D getGeometry() const;
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

    /// @name Cell
    /// @{
    /// @result The cell index containing the source. 
    /// @throws std::runtime_error if \c haveLocationInX() or
    ///         \c haveLocationInZ(), is false.
    [[nodiscard]] int getCell() const;
    /// @}

    /// @name Destructor
    /// @{
    /// @brief Destructor.
    ~Source2D();
    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @}
private:
    class Source2DImpl;
    std::unique_ptr<Source2DImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const Source2D &source);
}
#endif
