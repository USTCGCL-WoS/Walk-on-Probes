#pragma once

#include "solver/probe_set.hpp"
#include "solver/solver.hpp"
#include "solver/solver_api.h"

#include "core/math_defs.hpp"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
class WoPSolver : public Solver<ScalarType, DIM>
{
  public:
    using Solver<ScalarType, DIM>::Solver;

    ScalarType solve(const Vector<DIM>& targetPoint) override;

    void solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results) override;

    void configure(const nlohmann::json& j) override;

    void setMaxSteps(int n) { maxSteps_ = n; }
    void setWalksPerPixel(int n) { wpp_ = n; }
    void setEpsilon(double eps) { epsilon_ = eps; }
    void setRMin(double r) { rMin_ = r; }
    void setAlphaWalk(double a) { alphaWalk_ = a; }
    void setWMin(double w) { wMin_ = w; }

    void ensureProbes(const std::vector<Vector<DIM>>& targetPoints);

    const ProbeSet<ScalarType, DIM>& probeSet() const { return probes_; }

  private:
    ProbeSet<ScalarType, DIM> probes_;
    int maxSteps_ = 1024;
    int wpp_ = 1;
    double epsilon_ = EPSILON;
    double rMin_ = EPSILON;
    double alphaWalk_ = 0.7;
    double wMin_ = 1.0;
};

extern template class WoPSolver<Scalar<1>, 2>;
extern template class WoPSolver<Scalar<1>, 3>;
extern template class WoPSolver<Scalar<3>, 2>;
extern template class WoPSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
