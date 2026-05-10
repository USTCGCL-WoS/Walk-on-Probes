#include "solver/probe_set.hpp"
#include "kd_tree.hpp"
#include "probe_tree.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "pcg_random.hpp"
#include <spdlog/spdlog.h>

WOS_NAMESPACE_OPEN_SCOPE

namespace {

// Smoothstep covering weight: 3s^2 - 2s^3, where s = 1 - dist/radius
double
coveringWeight(double dist, double radius)
{
    if (dist >= radius)
        return 0.0;
    double s = 1.0 - dist / radius;
    return 3.0 * s * s - 2.0 * s * s * s;
}

template<int DIM>
struct Candidate
{
    Vector<DIM> point;
    double cumulativeWeight = 0.0;
    int origIdx = 0;
};

} // anonymous namespace

template<typename ScalarType, int DIM>
ProbeSet<ScalarType, DIM>::~ProbeSet()
{
    delete tree_;
}

// ==== build ====

template<typename ScalarType, int DIM>
void
ProbeSet<ScalarType, DIM>::build(const PoissonScene<ScalarType, DIM>& scene,
                                 const std::vector<Vector<DIM>>& targetPoints,
                                 uint64_t seed,
                                 double epsilon,
                                 double alphaRec,
                                 double alphaWalk,
                                 double wMin)
{
    alphaRec_ = alphaRec;
    alphaWalk_ = alphaWalk;
    wMin_ = wMin;

    probes_.clear();
    delete tree_;
    tree_ = nullptr;

    if (targetPoints.empty())
        return;

    // ==== Phase 1: Setup candidates ====

    std::vector<Candidate<DIM>> candidates;
    candidates.reserve(targetPoints.size());
    for (int i = 0; i < static_cast<int>(targetPoints.size()); ++i)
        candidates.push_back(Candidate<DIM>{ targetPoints[i], 0.0, i });

    // Shuffle
    pcg64 rng(seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);

    // Build KD-tree over candidate points for radius queries
    std::vector<Vector<DIM>> candidatePoints;
    candidatePoints.reserve(candidates.size());
    for (const auto& c : candidates)
        candidatePoints.push_back(c.point);

    using IndexType = size_t;
    using Adaptor = PointCloudAdaptor<DIM>;
    using KDTree =
      nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Adaptor>, Adaptor, DIM, IndexType>;

    Adaptor adaptor(candidatePoints);
    KDTree kdTree(DIM, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    kdTree.buildIndex();

    // ==== Phase 2: Greedy placement ====

    double alphaCov = std::min(alphaRec_, alphaWalk_);

    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        if (candidates[i].cumulativeWeight >= wMin_)
            continue;

        const auto& p = candidates[i].point;

        // Check all-boundary distance for probe radius
        fcpw::Interaction<DIM> boundaryInteract;
        if (!scene.findClosestBoundary(p, boundaryInteract))
            continue; // no boundary found, skip

        double radius = static_cast<double>(boundaryInteract.d) *
                        0.99; // Slightly shrink radius to avoid precision issues on boundary

        // Place probe
        {
            Probe<ScalarType, DIM> probe;
            probe.center = p;
            probe.radius = radius;
            probes_.push_back(std::move(probe));
        }

        // Update cumulative weights of affected remaining candidates
        double covRadius = alphaCov * radius;
        double covRadiusSq = covRadius * covRadius;

        std::vector<nanoflann::ResultItem<IndexType, double>> results;
        kdTree.radiusSearch(p.data(), covRadiusSq, results);

        std::vector<int> assigned; // For recording point indices assigned to this probe
        assigned.reserve(results.size());

        for (const auto& r : results) {
            int j = static_cast<int>(r.first);
            if (j > i) {
                // r.second is squared distance for L2 metric
                double dist = std::sqrt(r.second);
                candidates[j].cumulativeWeight += coveringWeight(dist, covRadius);
            }
            // If alphaRec_ == alphaCov, no need to search again
            if (alphaRec_ == alphaCov) {
                assigned.push_back(candidates[j].origIdx);
            }
        }

        // If alphaRec_ > alphaCov, search again using alphaRec_
        if (alphaRec_ > alphaCov) {
            double recRadius = alphaRec_ * radius;
            double recRadiusSq = recRadius * recRadius;

            results.clear();
            kdTree.radiusSearch(p.data(), recRadiusSq, results);
            assigned.reserve(results.size());

            for (const auto& r : results) {
                int j = static_cast<int>(r.first);
                assigned.push_back(candidates[j].origIdx);
            }
        }

        probes_.back().targetIndices = std::move(assigned);
    }

    // ==== Phase 3: Build spatial index ====

    if (!probes_.empty()) {
        // Compute bbox of probe effective regions
        Vector<DIM> bboxMin = probes_[0].center;
        Vector<DIM> bboxMax = probes_[0].center;
        for (const auto& probe : probes_) {
            double r = std::max(alphaRec_, alphaWalk_) * probe.radius;
            for (int d = 0; d < DIM; ++d) {
                bboxMin[d] = std::min(bboxMin[d], probe.center[d] - r);
                bboxMax[d] = std::max(bboxMax[d], probe.center[d] + r);
            }
        }

        tree_ = new detail::ProbeTree<ScalarType, DIM>(alphaRec_, alphaWalk_);
        tree_->build(probes_, bboxMin, bboxMax);
    }
}

// ==== findContainingProbe ====

template<typename ScalarType, int DIM>
const std::vector<const Probe<ScalarType, DIM>*>
ProbeSet<ScalarType, DIM>::findContainingProbes(const Vector<DIM>& p, double alpha) const
{
    if (!tree_ || probes_.empty())
        return {};

    return tree_->findContaining(probes_, p, alpha);
}

// ==== findOptimalProbe ====

template<typename ScalarType, int DIM>
const Probe<ScalarType, DIM>*
ProbeSet<ScalarType, DIM>::findOptimalProbe(const Vector<DIM>& p, double alpha) const
{
    if (!tree_ || probes_.empty())
        return nullptr;

    return tree_->findOptimal(probes_, p, alpha);
}

// ==== evaluate ====

template<typename ScalarType, int DIM>
const ScalarType
ProbeSet<ScalarType, DIM>::evaluate(const Vector<DIM>& p, int pointIdx, double kappa) const
{
    auto containing = findContainingProbes(p, alphaRec_);

    if (containing.empty()) {
        return ScalarType::NaN();
    }

    ScalarType result(0.0);
    double totalWeight = 0.0;

    for (const auto* probe : containing) {
        Vector<DIM> offset = p - probe->center;
        double d = offset.norm();
        double w = coveringWeight(d, probe->radius);
        if (w <= 0.0)
            continue;

        result += w * probe->evaluate(p, kappa);
        totalWeight += w;

        // Add cached source if pointIdx matches
        if (pointIdx >= 0) {
            for (size_t idx = 0; idx < probe->targetIndices.size(); ++idx) {
                if (probe->targetIndices[idx] == pointIdx) {
                    result += probe->srcAccum[idx];
                    break;
                }
            }
        }
    }

    return totalWeight > 0.0 ? result / totalWeight : ScalarType::NaN();
}

// ==== Explicit instantiations ====

template class WOS_SOLVER_API ProbeSet<Scalar<1>, 2>;
template class WOS_SOLVER_API ProbeSet<Scalar<1>, 3>;
template class WOS_SOLVER_API ProbeSet<Scalar<3>, 2>;
template class WOS_SOLVER_API ProbeSet<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
