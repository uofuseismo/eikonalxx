#ifndef EIKONALXX_RAY_PATH_3D_HPP
#define EIKONALXX_RAY_PATH_3D_HPP
#include <memory>
#include <vector>
namespace EikonalXX::Ray
{
 class Segment3D;
}
namespace EikonalXX::Ray
{
class Path3D
{
private:
    using PathType = std::vector<Segment3D>;
public:
    using iterator = typename PathType::iterator;
    using const_iterator = typename PathType::const_iterator;
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Path3D();
    /// @brief Copy constructor.
    /// @param[in] path  The path from which to initialize this class.
    Path3D(const Path3D &path);
    /// @brief Move constructor.
    /// @param[in,out] path  The path from which to initialize this class.
    ///                      On exit, path's behavior is undefined.
    Path3D(Path3D &&path) noexcept;
    /// @}

    /// @name Option 1 For Path Creation
    /// @{

    /// @brief Opens the path construction process.
    void open();
    /// @result True indicates the path is open for construction.
    [[nodiscard]] bool isOpen() const noexcept;
    /// @brief Appends to the path.
    /// @param[in] segment  The segment to append to the path.
    /// @throws std::invalid_argument if the segment does not have a start/end
    ///         point and a velocity.  Additionally, if this is the not the
    ///         source segment then the segment's start point must coincide
    ///         with the previous segment's end point.
    /// @throws std::runtime_error if \c isOpen() is false.
    void append(const Segment3D &segment);
    /// @brief Appends to the path.
    /// @param[in,out] segment  The segment to append to the path.
    ///                         On exit, segment's behavior is undefined.
    void append(Segment3D &&segment);
    /// @brief Closes the path construction process.
    /// @throws std::runtime_error if \c isOpen() is false.
    void close();
    /// @}

    /// @name Option 2 For Path Creation
    /// @{

    /// @brief Sets the ray path.
    /// @param[in] segments  The segments to set.
    /// @throws std::invalid_argument if the i'th segment's end point does
    ///         not equal the i+1'th segment's start point, or the velocity
    ///         is not set on a segment.
    void set(const std::vector<Segment3D> &segments);
    /// @brief Sets the ray path.
    /// @param[in,out] segments  The segments to set.  On exit, segments
    ///                          behavior is undefined.
    /// @throws std::invalid_argument if the i'th segment's end point does
    ///         not equal the i+1'th segment's start point, or the velocity
    ///         is not set on a segment.
    void set(std::vector<Segment3D> &&segments);
    /// @}

    /// @result The number of segments. 
    [[nodiscard]] size_t size() const noexcept;
    /// @result The total travel time in seconds along the ray.
    [[nodiscard]] double getTravelTime() const;
    /// @result The total length of the ray in meters.
    [[nodiscard]] double getLength() const;

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Path3D();
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] path   The path to copy to this.
    /// @result A deep copy of the path. 
    Path3D &operator=(const Path3D &path);
    /// @brief Move assignment operator.
    /// @result The memory from path moved to this.
    Path3D &operator=(Path3D &&path) noexcept;
    /// @result A reference to the first ray segment in the path.
    iterator begin();
    /// @result A constant reference to the first ray segment in the path.
    constexpr const_iterator begin() const noexcept;
    /// @result A constant reference to the first ray segment in the path.
    constexpr const_iterator cbegin() const noexcept;

    /// @result A reference to the last ray segment in the path. 
    iterator end(); 
    /// @result A reference to the last ray segment in the path.
    constexpr const_iterator end() const noexcept;
    /// @result A reference to the last ray segment in the path.
    constexpr const_iterator cend() const noexcept;

    /// @param[in] index  The index of the desired segment.
    /// @result A reference to the segment at the given position.
    /// @throws std::invalid_argument if this is out of bounds.
    Segment3D& at(size_t index);
    /// @param[in] index  The index of the desired segment.
    /// @result A reference to the segment at the given position.
    /// @throws std::invalid_argument if this is out of bounds.
    const Segment3D& at(size_t index) const;
    /// @}
private:
    class Path3DImpl;
    std::unique_ptr<Path3DImpl> pImpl;
};
}
#endif
