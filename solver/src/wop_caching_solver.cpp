#include "solver/wop_caching_solver.hpp"
#include "core/green.hpp"
#include "core/math_defs.hpp"
#include "core/sampling_uniform.hpp"
#include "core/vector.hpp"
#include "fcpw/core/interaction.h"
#include "solver/walk_path.hpp"
#include "walk.hpp"
#include <cstddef>
#include <nlohmann/json.hpp>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

WOS_NAMESPACE_OPEN_SCOPE

// ==== configure ====

template<typename ScalarType, int DIM>
void
WoPCachingSolver<ScalarType, DIM>::configure(const nlohmann::json& j)
{
    if (j.contains("max_steps"))
        setMaxSteps(j["max_steps"].get<int>());
    if (j.contains("epsilon"))
        setEpsilon(j["epsilon"].get<double>());
    if (j.contains("r_min"))
        setRMin(j["r_min"].get<double>());
    if (j.contains("num_iterations"))
        setNumIterations(j["num_iterations"].get<int>());
    if (j.contains("fourier_modes"))
        setFourierModes(j["fourier_modes"].get<int>());
    if (j.contains("sh_modes"))
        setSHModes(j["sh_modes"].get<int>());
    if (j.contains("lambda"))
        setBoundarySampleFactor(j["lambda"].get<double>());
    if (j.contains("mu"))
        setSourceSampleFactor(j["mu"].get<double>());
    if (j.contains("min_boundary_samples"))
        setMinBoundarySamples(j["min_boundary_samples"].get<int>());
    if (j.contains("min_source_samples"))
        setMinSourceSamples(j["min_source_samples"].get<int>());
    if (j.contains("alpha_rec"))
        setAlphaRec(j["alpha_rec"].get<double>());
    if (j.contains("alpha_walk"))
        setAlphaWalk(j["alpha_walk"].get<double>());
    if (j.contains("w_min"))
        setWMin(j["w_min"].get<double>());
    if (j.contains("wost_wpp"))
        setWoStWalksPerPixel(j["wost_wpp"].get<int>());
    if (j.contains("use_sn"))
        setUseSN(j["use_sn"].get<bool>());
    if (j.contains("use_cv"))
        setUseCV(j["use_cv"].get<bool>());
}

// ==== runIteration ====
template<typename ScalarType, int DIM>
void
WoPCachingSolver<ScalarType, DIM>::runIteration(int iter)
{
    auto& probesVec = probes_.probes();
    const int numProbes = static_cast<int>(probesVec.size());
    if (numProbes == 0)
        return;

    const auto& scene = this->scene();
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;

    // ---- Phase 1: sample boundaries and regenerate starting points for WoP ----
    if (iter == 0) {
        allBoundaryPoints_.clear();
        std::vector<Vector<DIM>> boundarySamples;
        for (int p = 0; p < numProbes; ++p) {
            auto& probe = probesVec[p];
            int N = std::max(static_cast<int>(lambda_ * std::pow(probe.radius, DIM - 1)), minBoundarySamples_);
            sampleBoundarySamples(probe, N, boundarySamples);
            allBoundaryPoints_.insert(allBoundaryPoints_.end(), boundarySamples.begin(), boundarySamples.end());
            boundaryPointToProbeIdx_.insert(boundaryPointToProbeIdx_.end(), boundarySamples.size(), p);
        }
    }

    // ---- Phase 2: WoP walks + reverse accumulation ----
    if (!allBoundaryPoints_.empty()) {
        spdlog::info("WoPCachingSolver: running WoP walks for {} boundary samples (iteration {})...",
                     allBoundaryPoints_.size(),
                     iter);
#pragma omp parallel
        {
            WalkPath<ScalarType, DIM> localPath;
            fcpw::Interaction<static_cast<size_t>(DIM)> interaction;

#pragma omp for schedule(dynamic)
            for (size_t k = 0; k < allBoundaryPoints_.size(); ++k) {
                localPath.clear();
                walkWoP(scene,
                        probes_,
                        localPath,
                        allBoundaryPoints_[k],
                        maxSteps_,
                        epsilon_,
                        rMin_,
                        this->rng_,
                        this->dist01_,
                        interaction);
                if (localPath.dirichlet.isNaN()) {
                    --k; // retry this walk
                    continue;
                }

                // Reverse accumulation onto probes
                const auto& pathPositions = localPath.positions;
                const auto& pathSourceContribs = localPath.sourceContribs;
                const auto& pathNeumannContribs = localPath.neumannContribs;
                const auto& pathPdfs = localPath.pdfs;
                const auto& pathProbeIndices = localPath.probeIndices;
                const int pathCount = localPath.count;

                ScalarType estimate = localPath.dirichlet;
                Vector<DIM> pos = localPath.dirichletPos;
                double pdf;
                int probeIdx;
                if (pathCount > 0) {
                    for (int j = pathCount - 1; j >= 0; --j) {
                        if (j != pathCount - 1) {
                            estimate += pathSourceContribs[j + 1] + pathNeumannContribs[j + 1];
                            pos = pathPositions[j + 1];
                        }
                        pdf = pathPdfs[j];
                        probeIdx = pathProbeIndices[j];
                        if (probeIdx >= 0) {
                            auto& probe = probesVec[probeIdx];
                            probe.addSampleUnlocked(pos, estimate, pdf);
                        }
                    }
                    // Add contribution to probe of starting point
                    estimate += pathSourceContribs[0] + pathNeumannContribs[0];
                    pos = pathPositions[0];
                }
                pdf = (DIM == 2) ? 0.5 * INV_PI : 0.25 * INV_PI;
                probeIdx = boundaryPointToProbeIdx_[k];
                auto& probe = probesVec[probeIdx];
                probe.addSampleUnlocked(pos, estimate, pdf);
            }
        }
    }

    // ---- Phase 4: source sampling per probe ----
    spdlog::info("WoPCachingSolver: sampling sources for {} probes (iteration {})...", numProbes, iter);
#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < numProbes; ++p) {
        auto& probe = probesVec[p];
        if (probe.targetIndices.empty()) {
            probe.srcCount++;
            continue;
        }

        double R = probe.radius;

        int M = std::max(static_cast<int>(mu_ * std::pow(R, DIM)), minSourceSamples_);

        for (int j = 0; j < M; ++j) {
            auto ySample = sampleUniformBall<DIM>(this->rng_, this->dist01_, R);
            Vector<DIM> Yj = probe.center + ySample.point;
            ScalarType fVal = scene.source(Yj);

            for (size_t localIdx = 0; localIdx < probe.targetIndices.size(); ++localIdx) {
                int t = probe.targetIndices[localIdx];
                Vector<DIM> xRel = targetPoints_[t] - probe.center;
                double G = greensFunction<DIM>(xRel, ySample.point, R, kappa);
                probe.srcAccum[localIdx] += G * fVal / (ySample.pdf * M);
            }
        }

        probe.srcCount++;
    }
}

// ==== sampleBoundarySamples (internal helper) ====
template<typename ScalarType, int DIM>
void
WoPCachingSolver<ScalarType, DIM>::sampleBoundarySamples(const Probe<ScalarType, DIM>& probe,
                                                         int N,
                                                         std::vector<Vector<DIM>>& boundarySamples)
{
    const auto& scene = this->scene();
    boundarySamples.clear();
    boundarySamples.reserve(N);

    if constexpr (DIM == 2) {
        for (int i = 0; i < N; ++i) {
            double u = this->dist01_(this->rng_);
            double theta = 2.0 * PI * (static_cast<double>(i) + u) / static_cast<double>(N);

            // Boundary sample position (absolute)
            Vector<DIM> dir;
            dir << std::cos(theta), std::sin(theta);
            Vector<DIM> boundaryPoint = probe.center + probe.radius * dir;
            if (scene.isInteriorOnly() && !scene.isInsideObject(boundaryPoint)) {
                // Regenerate direction until we get a valid boundary point
                --i;
                continue;
            }
            boundarySamples.push_back(boundaryPoint);
        }
    } else if constexpr (DIM == 3) {
        // Equal-area stratification on sphere: factor N into nTheta * nPhi
        int nPhi = std::max(1, static_cast<int>(std::round(std::sqrt(static_cast<double>(N)))));
        int nTheta = (N + nPhi - 1) / nPhi;
        int sampleIdx = 0;

        for (int ti = 0; ti < nTheta && sampleIdx < N; ++ti) {
            double uT = this->dist01_(this->rng_);
            double cosColat = 1.0 - 2.0 * (static_cast<double>(ti) + uT) / static_cast<double>(nTheta);
            cosColat = std::clamp(cosColat, -1.0, 1.0);
            double sinColat = std::sqrt(std::max(0.0, 1.0 - cosColat * cosColat));

            for (int pi = 0; pi < nPhi && sampleIdx < N; ++pi) {
                double uP = this->dist01_(this->rng_);
                double phi = 2.0 * PI * (static_cast<double>(pi) + uP) / static_cast<double>(nPhi);

                // Boundary sample position
                Vector<DIM> dir;
                dir << sinColat * std::cos(phi), sinColat * std::sin(phi), cosColat;
                Vector<DIM> boundaryPoint = probe.center + probe.radius * dir;
                if (scene.isInteriorOnly() && !scene.isInsideObject(boundaryPoint)) {
                    // Regenerate this sample until we get a valid boundary point
                    --pi;
                    continue;
                }

                ++sampleIdx;
                boundarySamples.push_back(boundaryPoint);
            }
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
}

// ==== buildCache ====

template<typename ScalarType, int DIM>
void
WoPCachingSolver<ScalarType, DIM>::buildCache()
{
    const auto& scene = this->scene();

    if (targetPoints_.empty()) {
        spdlog::warn("WoPCachingSolver::buildCache: no target points set, skipping");
        return;
    }

    // 1. Build probe set
    probes_.build(scene, targetPoints_, this->seed_, epsilon_, alphaRec_, alphaRec_, wMin_);

    if (probes_.empty()) {
        spdlog::warn("WoPCachingSolver::buildCache: no probes placed, skipping");
        return;
    }

    spdlog::info(
      "WoPCachingSolver: {} probes placed for {} target points", probes_.probes().size(), targetPoints_.size());

    // 2. Initialize expansion structures
    for (auto& probe : probes_.probes()) {
        if constexpr (DIM == 2) {
            probe.initExpansion(fourierModes_);
        } else if constexpr (DIM == 3) {
            probe.initExpansion(shModes_);
        } else {
            static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        }
        probe.srcAccum.assign(probe.targetIndices.size(), ScalarType(0.0));
    }

    // 3. Run iterations
    spdlog::info("WoPCachingSolver: running {} iterations...", numIterations_);
    for (int iter = 0; iter < numIterations_; ++iter) {
        runIteration(iter);
        if ((iter + 1) % 10 == 0 || iter == numIterations_ - 1)
            spdlog::info("  iteration {}/{}", iter + 1, numIterations_);
    }

    // 4. Finalize coefficients
    for (auto& probe : probes_.probes())
        probe.finalizeCoefficients(useSN_, useCV_);

    cacheBuilt_ = true;
    spdlog::info("WoPCachingSolver: cache built successfully");
}

// ==== solvePoint (internal reconstruction helper) ====

template<typename ScalarType, int DIM>
ScalarType
WoPCachingSolver<ScalarType, DIM>::solvePoint(const Vector<DIM>& p,
                                              int pointIdx,
                                              pcg64& rng,
                                              std::uniform_real_distribution<double>& dist,
                                              fcpw::Interaction<static_cast<size_t>(DIM)>& interaction)
{
    const auto& scene = this->scene();
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;

    ScalarType result = probes_.evaluate(p, pointIdx, kappa);
    if (result.isNaN()) {
        // Fallback: averaged WoSt
        ScalarType sum(0.0);
        for (int w = 0; w < wostWpp_; ++w) {
            auto walkResult = walkWoSt(scene, p, maxSteps_, epsilon_, rMin_, rng, dist, interaction);
            if (walkResult.isNaN()) {
                --w; // retry this walk
                continue;
            }
            sum += walkResult;
        }
        return sum / static_cast<double>(wostWpp_);
    }
    return result;
}

// ==== solve (single point) ====

template<typename ScalarType, int DIM>
ScalarType
WoPCachingSolver<ScalarType, DIM>::solve(const Vector<DIM>& targetPoint)
{
    fcpw::Interaction<static_cast<size_t>(DIM)> interaction;
    return solvePoint(targetPoint, -1, this->rng_, this->dist01_, interaction);
}

// ==== solve (batch) ====

template<typename ScalarType, int DIM>
void
WoPCachingSolver<ScalarType, DIM>::solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results)
{
    if (!cacheBuilt_) {
        targetPoints_ = points;
        buildCache();
    }

    results.resize(points.size());

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto threadRng = Solver<ScalarType, DIM>::makeRngFromSeed(this->seed_, 1 + tid);
        std::uniform_real_distribution<double> threadDist(0.0, 1.0);
        fcpw::Interaction<static_cast<size_t>(DIM)> interaction;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < static_cast<int>(points.size()); ++i) {
            results[i] = solvePoint(points[i], i, threadRng, threadDist, interaction);
        }
    }
}

template class WOS_SOLVER_API WoPCachingSolver<Scalar<1>, 2>;
template class WOS_SOLVER_API WoPCachingSolver<Scalar<1>, 3>;
template class WOS_SOLVER_API WoPCachingSolver<Scalar<3>, 2>;
template class WOS_SOLVER_API WoPCachingSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
