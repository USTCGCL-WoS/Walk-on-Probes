#pragma once

#include <cstdint>
#include <vector>

#include "solver/probe.hpp"
#include "solver/solver_api.h"

#include "core/scalar.hpp"
#include "core/vector.hpp"
#include "scene/poisson.hpp"

WOS_NAMESPACE_OPEN_SCOPE

namespace detail {
template<typename, int>
struct ProbeTree;
} // namespace detail

template<typename ScalarType, int DIM>
class ProbeSet
{
  public:
    void build(const PoissonScene<ScalarType, DIM>& scene,
               const std::vector<Vector<DIM>>& targetPoints,
               uint64_t seed,
               double epsilon,
               double alphaRec,
               double alphaWalk,
               double wMin);

    const std::vector<const Probe<ScalarType, DIM>*> findContainingProbes(const Vector<DIM>& p, double alpha) const;

    const Probe<ScalarType, DIM>* findOptimalProbe(const Vector<DIM>& p, double alpha) const;

    const ScalarType evaluate(const Vector<DIM>& p, int pointIdx, double kappa) const;

    ~ProbeSet();

    std::vector<Probe<ScalarType, DIM>>& probes() { return probes_; }
    const std::vector<Probe<ScalarType, DIM>>& probes() const { return probes_; }
    bool empty() const { return probes_.empty(); }

  private:
    std::vector<Probe<ScalarType, DIM>> probes_;
    detail::ProbeTree<ScalarType, DIM>* tree_ =
      nullptr; // raw ptr: destructor defined in .cpp where ProbeTree is complete

    double alphaRec_ = 0.9;
    double alphaWalk_ = 0.7;
    double wMin_ = 1.0;
};

extern template class ProbeSet<Scalar<1>, 2>;
extern template class ProbeSet<Scalar<1>, 3>;
extern template class ProbeSet<Scalar<3>, 2>;
extern template class ProbeSet<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
