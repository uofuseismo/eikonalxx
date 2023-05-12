#ifndef EIKONALXX_STATION_2D_HPP
#define EIKONALXX_STATION_2D_HPP
#include <ostream>
#include <memory>
namespace EikonalXX
{
class Geometry2D;
/// @class Station2D "station2d.hpp" "eikonalxx/station2d.hpp"
/// @brief Defines a station in a 2D model.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Station2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Station2D();
    /// @brief Copy constructor.
    /// @param[in] station  The station location from which to initialize
    ///                     this class.
    Station2D(const Station2D &station);
    /// @brief Move constructor.
    /// @param[in,out] station  The station location from which to initialize
    ///                         this class.  On exit, station's behavior is
    ///                         undefined.
    Station2D(Station2D &&station) noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] station  The station location to copy to this.
    /// @result A deep copy of the input station.
    Station2D& operator=(const Station2D &station);
    /// @brief Move assignment operator.
    /// @param[in] station  The station location whose memory will be moved to
    ///                     this.  On exit, station's behavior is undefined.
    /// @result The memory from station moved to this.
    Station2D& operator=(Station2D &&station) noexcept;
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

    /// @name Name
    /// @{

    /// @brief Sets the station's name.
    /// @param[in] name  The name of the station.
    void setName(const std::string &name) noexcept;
    /// @result The station's name.
    std::string getName() const noexcept;
    /// @}

    /// @name X Location
    /// @{

    /// @param[in] x   The station location in x in meters.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    /// @throws std::invalid_argument if the x position is not in the model.
    void setLocationInX(double x);
    /// @result The station's location in x in meters.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] double getLocationInX() const;
    /// @result The station's offset from the origin in x in meters.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] double getOffsetInX() const;
    /// @result The x cell containing the station.
    /// @throws std::runtime_error if \c haveLocationInX() is false.
    [[nodiscard]] int getCellInX() const;
    /// @result True indicates that the x location was set.
    [[nodiscard]] bool haveLocationInX() const noexcept;
    /// @}

    /// @name Z Location
    /// @{

    /// @param[in] z    The station location in z in meters.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    /// @throws std::invalid_argument if the z position is not in the model.
    void setLocationInZ(double z);
    /// @brief Sets the z position to the free surface in the rectilinear
    ///        geometry.  This is useful when computing reciprocity travel time
    ///        fields or modeling surface shots.
    /// @throws std::runtime_error if \c haveGeometry() is false.
    void setZToFreeSurface();
    /// @result The station's location in z in meters.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] double getLocationInZ() const;
    /// @result The station's offset from the origin in z in meters.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] double getOffsetInZ() const;
    /// @result The z cell containing the station.
    /// @throws std::runtime_error if \c haveLocationInZ() is false.
    [[nodiscard]] int getCellInZ() const;
    /// @result True indicates that the z location was set.
    [[nodiscard]] bool haveLocationInZ() const noexcept;
    /// @}

    /// @name Cell
    /// @{

    /// @result The cell index containing the station. 
    /// @throws std::runtime_error if \c haveLocationInX() or
    ///         \c haveLocationInZ(), is false.
    [[nodiscard]] int getCell() const;
    /// @}

    /// @name Destructor
    /// @{

    /// @brief Destructor.
    ~Station2D();
    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @}
private:
    class Station2DImpl;
    std::unique_ptr<Station2DImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const Station2D &station);
}
#endif
