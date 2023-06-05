#ifndef PYEIKONALXX_SOLVER_2D_HPP
#define PYEIKONALXX_SOLVER_2D_HPP
#include <memory>
#include <vector>
#include <pybind11/pybind11.h>
namespace EikonalXX
{
 template<class T> class Solver2D;
}
namespace PEikonalXX
{
 class Geometry2D;
 class Model2D;
 class SolverOptions;
 class Source2D;
 class Station2D;
}
namespace PEikonalXX
{
class Solver2D
{
public:
    Solver2D();

    void initialize(const SolverOptions &options,
                    const Geometry2D &geometry);
    [[nodiscard]] SolverOptions getOptions() const;
    [[nodiscard]] Geometry2D getGeometry() const;
    [[nodiscard]] bool isInitialized() const noexcept;

    void setVelocityModel(const Model2D &velocityModel);
    [[nodiscard]] Model2D getVelocityModel() const;
    [[nodiscard]] bool haveVelocityModel() const noexcept;
    [[nodiscard]] double getSlowness(const int iCellX, const int iCellZ) const;

    void setSource(const Source2D &source);
    void setSource(const std::pair<double, double> &location);
    [[nodiscard]] Source2D getSource() const;
    [[nodiscard]] bool haveSource() const noexcept;

    void solve(); 

    [[nodiscard]] std::vector<double> getTravelTimeField() const;
    [[nodiscard]] bool haveTravelTimeField() const noexcept;
    [[nodiscard]] std::vector<double> getTravelTimeGradientField() const;
    [[nodiscard]] bool haveTravelTimeGradientField() const noexcept;

    ~Solver2D();
private:
    std::unique_ptr<EikonalXX::Solver2D<double>> mSolver;
};
void initializeSolver2D(pybind11::module &m);
}
#endif
