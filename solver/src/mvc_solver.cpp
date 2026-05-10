#include "solver/mvc_solver.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "core/green.hpp"
#include "core/math_defs.hpp"
#include "core/sampling_uniform.hpp"
#include "core/vector.hpp"
#include "fcpw/core/interaction.h"
#include "scene/poisson.hpp"
#include "walk.hpp"

#include <nlohmann/json.hpp>

#include <spdlog/spdlog.h>

#include <omp.h>

#include "kd_tree.hpp"

#include <nanoflann.hpp>

WOS_NAMESPACE_OPEN_SCOPE

namespace {

template<int DIM>
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor<DIM>>,
                                                   PointCloudAdaptor<DIM>,
                                                   DIM,
                                                   size_t>;

template<typename ScalarType, int DIM>
ScalarType
queryOnePoint(const PoissonScene<ScalarType, DIM>& scene,
              const Vector<DIM>& p,
              double epsilon,
              int wostWpp,
              int maxSteps,
              double rMin,
              int minSourceSamples,
              double sourceRatio,
              pcg64& rng,
              std::uniform_real_distribution<double>& dist01,
              fcpw::Interaction<static_cast<size_t>(DIM)>& interaction,
              const std::vector<Vector<DIM>>& cachePositions,
              const std::vector<ScalarType>& cacheSolutions,
              const std::vector<double>& cacheWeights,
              KDTree<DIM>& kdTree)
{
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;

    // Find distance to closest Dirichlet boundary
    double distD = std::numeric_limits<double>::max();
    if (scene.findClosestDirichletBoundary(p, interaction))
        distD = static_cast<double>(interaction.d);

    // Directly on Dirichlet boundary → return BC value
    if (distD < epsilon) {
        auto bc = scene.boundaryCondition(interaction, p, BoundaryType::Dirichlet);
        return bc.a;
    }

    // Distance to interest boundary
    double distInterest = scene.distanceToInterestBoundary(p, interaction);
    if (distInterest < 0.0)
        return ScalarType::NaN();

    // Ball radius
    double R = std::min(distD, distInterest) * 0.99;

    // Range-search KD-tree for cache points within ball
    double radiusSq = R * R;
    std::vector<nanoflann::ResultItem<size_t, double>> results;
    kdTree.radiusSearch(p.data(), radiusSq, results);

    // Fallback: no cache points in ball → walkWoSt
    if (results.empty()) {
        ScalarType sum(0.0);
        for (int w = 0; w < wostWpp; ++w) {
            auto walkResult = walkWoSt(scene, p, maxSteps, epsilon, rMin, rng, dist01, interaction);
            if (walkResult.isNaN()) {
                --w;
                continue;
            }
            sum += walkResult;
        }
        return sum / static_cast<double>(wostWpp);
    }

    // Weighted average of cache point solutions
    ScalarType numerator(0.0);
    double weightSum = 0.0;

    for (const auto& result : results) {
        size_t idx = result.first;
        Vector<DIM> offset = cachePositions[idx] - p;
        double normConst = 1.0f;
        if (kappa > EPSILON) {
            double kR = kappa * R;
            if constexpr (DIM == 2) {
                normConst = 0.5 * kR / std::cyl_bessel_i(1, kR);
            } else if constexpr (DIM == 3) {
                normConst = kR * kR * kR / (3.0 * (kR * std::cosh(kR) - std::sinh(kR)));
            } else {
                static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
            }
        }

        double w = cacheWeights[idx];

        numerator += w * cacheSolutions[idx] * normConst;
        weightSum += w;
    }

    ScalarType boundaryPart = (weightSum > 0.0) ? numerator / weightSum : ScalarType::NaN();
    if (boundaryPart.isNaN()) {
        ScalarType sum(0.0);
        for (int w = 0; w < wostWpp; ++w) {
            auto walkResult = walkWoSt(scene, p, maxSteps, epsilon, rMin, rng, dist01, interaction);
            if (walkResult.isNaN()) {
                --w;
                continue;
            }
            sum += walkResult;
        }
        return sum / static_cast<double>(wostWpp);
    }

    // Source term MC integration
    int M = std::max(minSourceSamples, static_cast<int>(static_cast<double>(results.size()) * sourceRatio));
    ScalarType sourceAccum(0.0);

    for (int j = 0; j < M; ++j) {
        auto sample = sampleUniformBall<DIM>(rng, dist01, R);
        Vector<DIM> y = sample.point;
        ScalarType fVal = scene.source(p + y);
        double gamma = greensFunctionAtCenter(y, R);
        if constexpr (DIM == 2) {
            gamma += 0.25 * (y.squaredNorm() - R * R) / (PI * R * R);
        } else if constexpr (DIM == 3) {
            gamma += 0.125 * (y.squaredNorm() - R * R) / (PI * R * R * R);
        } else {
            static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        }
        sourceAccum += gamma * fVal / sample.pdf;
    }
    sourceAccum /= static_cast<double>(M);

    return boundaryPart + sourceAccum;
}

} // namespace

// ==== configure ====

template<typename ScalarType, int DIM>
void
MVCSolver<ScalarType, DIM>::configure(const nlohmann::json& j)
{
    if (j.contains("num_iterations"))
        setNumIterations(j["num_iterations"].get<int>());
    if (j.contains("num_cache_points"))
        setNumCachePoints(j["num_cache_points"].get<int>());
    if (j.contains("wost_wpp"))
        setWoStWalksPerPixel(j["wost_wpp"].get<int>());
    if (j.contains("max_steps"))
        setMaxSteps(j["max_steps"].get<int>());
    if (j.contains("epsilon"))
        setEpsilon(j["epsilon"].get<double>());
    if (j.contains("r_min"))
        setRMin(j["r_min"].get<double>());
    if (j.contains("source_ratio"))
        setSourceRatio(j["source_ratio"].get<double>());
    if (j.contains("min_source_samples"))
        setMinSourceSamples(j["min_source_samples"].get<int>());
    if (j.contains("mu"))
        setMu(j["mu"].get<double>());
    if (j.contains("lambda"))
        setLambda(j["lambda"].get<double>());
    if (j.contains("adapt_q"))
        setAdaptQ(j["adapt_q"].get<double>());
}

// ==== generateCachePoints ====

template<typename ScalarType, int DIM>
void
MVCSolver<ScalarType, DIM>::generateCachePoints(std::vector<CachePoint>& out)
{
    const auto& scene = this->scene();

    Vector<DIM> renderMin = scene.getRenderMin();
    Vector<DIM> renderMax = scene.getRenderMax();

    out.clear();
    out.reserve(numCachePoints_);

    fcpw::Interaction<static_cast<size_t>(DIM)> interaction;

    while (static_cast<int>(out.size()) < numCachePoints_) {
        // Sample uniformly in ROI bounding box
        Vector<DIM> p;
        for (int d = 0; d < DIM; ++d)
            p[d] = renderMin[d] + this->dist01_(this->rng_) * (renderMax[d] - renderMin[d]);

        // Check ROI membership
        if (!scene.isInsideROI(p))
            continue;

        double importanceWeight = 1.0;

        if (adaptQ_ > 0.0) {
            double acceptProb = 1.0 / adaptQ_;

            if (!scene.findClosestBoundary(p, interaction))
                continue;

            double r_d = static_cast<double>(interaction.d);
            double r_i = scene.distanceToInterestBoundary(p, interaction);
            double r = std::min(r_d, r_i);

            if (r > mu_) {
                double factor;
                if constexpr (DIM == 2) {
                    factor = mu_ * mu_ / (r * r);
                } else if constexpr (DIM == 3) {
                    factor = mu_ * mu_ * mu_ / (r * r * r);
                } else {
                    static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
                }
                acceptProb = std::max(factor / adaptQ_, lambda_ / adaptQ_);
            }

            if (this->dist01_(this->rng_) > acceptProb) {
                continue;
            }

            importanceWeight = 1.0 / acceptProb;
        }

        out.push_back({ p, ScalarType(0.0), importanceWeight });
    }

    spdlog::info("MVCSolver: generated {} cache points (target: {})", out.size(), numCachePoints_);
}

// ==== buildCache ====

template<typename ScalarType, int DIM>
void
MVCSolver<ScalarType, DIM>::buildCache()
{
    const auto& scene = this->scene();

    if (targetPoints_.empty()) {
        spdlog::warn("MVCSolver::buildCache: no target points set, skipping");
        return;
    }

    // Initialize running averages
    accumulatedSolutions_.assign(targetPoints_.size(), ScalarType(0.0));
    accumulatedCounts_.assign(targetPoints_.size(), 0);

    // Create delegate WoSt solver (reused across rounds)
    woStSolver_ = std::make_unique<WoStSolver<ScalarType, DIM>>(scene, this->seed_ + 1);
    woStSolver_->setMaxSteps(maxSteps_);
    woStSolver_->setEpsilon(epsilon_);
    woStSolver_->setRMin(rMin_);
    woStSolver_->setWalksPerPixel(wostWpp_);

    spdlog::info("MVCSolver: building cache for {} target points over {} rounds, {} cache points/round",
                 targetPoints_.size(),
                 numIterations_,
                 numCachePoints_);

    for (int round = 0; round < numIterations_; ++round) {
        // a. Generate cache points
        std::vector<CachePoint> cachePoints;
        generateCachePoints(cachePoints);
        if (cachePoints.empty()) {
            spdlog::warn("MVCSolver: no cache points generated in round {}, skipping", round);
            continue;
        }

        // b. Extract positions, solve with WoSt
        std::vector<Vector<DIM>> cachePositions;
        std::vector<double> cacheWeights;
        cachePositions.reserve(cachePoints.size());
        cacheWeights.reserve(cachePoints.size());
        for (const auto& cp : cachePoints) {
            cachePositions.push_back(cp.position);
            cacheWeights.push_back(cp.importanceWeight);
        }

        std::vector<ScalarType> cacheSolutions;
        woStSolver_->solve(cachePositions, cacheSolutions);

        // c-d. Build KD-tree
        PointCloudAdaptor<DIM> adaptor(cachePositions);
        KDTree<DIM> kdTree(DIM, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        kdTree.buildIndex();

        // e. Query each target point in parallel
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto threadRng = Solver<ScalarType, DIM>::makeRngFromSeed(this->seed_, 2 + round + tid);
            std::uniform_real_distribution<double> threadDist(0.0, 1.0);
            fcpw::Interaction<static_cast<size_t>(DIM)> interaction;

#pragma omp for schedule(dynamic)
            for (int i = 0; i < static_cast<int>(targetPoints_.size()); ++i) {
                ScalarType result = queryOnePoint<ScalarType, DIM>(scene,
                                                                   targetPoints_[i],
                                                                   epsilon_,
                                                                   wostWpp_,
                                                                   maxSteps_,
                                                                   rMin_,
                                                                   minSourceSamples_,
                                                                   sourceRatio_,
                                                                   threadRng,
                                                                   threadDist,
                                                                   interaction,
                                                                   cachePositions,
                                                                   cacheSolutions,
                                                                   cacheWeights,
                                                                   kdTree);

                if (!result.isNaN()) {
                    accumulatedSolutions_[i] += result;
                    accumulatedCounts_[i]++;
                }
            }
        }

        if ((round + 1) % 10 == 0 || round == numIterations_ - 1 || round == 0)
            spdlog::info("  round {}/{}", round + 1, numIterations_);
    }

    // Average across rounds
    for (size_t i = 0; i < targetPoints_.size(); ++i) {
        if (accumulatedCounts_[i] > 0)
            accumulatedSolutions_[i] /= static_cast<double>(accumulatedCounts_[i]);
        else
            accumulatedSolutions_[i] = ScalarType::NaN();
    }

    cacheBuilt_ = true;
    spdlog::info("MVCSolver: cache built successfully");
}

// ==== solve (single point) ====

template<typename ScalarType, int DIM>
ScalarType
MVCSolver<ScalarType, DIM>::solve(const Vector<DIM>& targetPoint)
{
    if (!cacheBuilt_) {
        targetPoints_ = { targetPoint };
        buildCache();
    }
    return accumulatedSolutions_[0];
}

// ==== solve (batch) ====

template<typename ScalarType, int DIM>
void
MVCSolver<ScalarType, DIM>::solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results)
{
    if (!cacheBuilt_) {
        targetPoints_ = points;
        buildCache();
    }
    results = accumulatedSolutions_;
}

// ==== Explicit instantiations ====

template class WOS_SOLVER_API MVCSolver<Scalar<1>, 2>;
template class WOS_SOLVER_API MVCSolver<Scalar<1>, 3>;
template class WOS_SOLVER_API MVCSolver<Scalar<3>, 2>;
template class WOS_SOLVER_API MVCSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
