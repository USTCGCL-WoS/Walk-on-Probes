#pragma once

#include "core/green.hpp"
#include "core/math_defs.hpp"
#include "core/poisson_kernel.hpp"
#include "core/sample.hpp"
#include "core/sampling_green.hpp"
#include "core/sampling_poisson_kernel.hpp"
#include "core/sampling_uniform.hpp"
#include "core/scalar.hpp"
#include "core/vector.hpp"
#include "scene/poisson.hpp"
#include "solver/probe_set.hpp"
#include "solver/walk_path.hpp"
#include <spdlog/spdlog.h>

#include "pcg_random.hpp"

#include <cassert>
#include <limits>
#include <random>

WOS_NAMESPACE_OPEN_SCOPE

// ==== Walk on Spheres ====
template<typename ScalarType, int DIM>
ScalarType
walkWoS(const PoissonScene<ScalarType, DIM>& scene,
        const Vector<DIM>& startPoint,
        int maxSteps,
        double epsilon,
        pcg64& rng,
        std::uniform_real_distribution<double>& dist01,
        fcpw::Interaction<static_cast<size_t>(DIM)>& interaction)
{
    Vector<DIM> pos = startPoint;
    ScalarType result(0.0);
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;
    double screeningThroughput = 1.0;

    for (int step = 0; step < maxSteps; ++step) {
        // Interior point checks
        if (scene.isInteriorOnly() && !scene.isInsideObject(pos)) {
            // Walked outside of interior domain, return NaN
            return ScalarType::NaN();
        }

        if (!scene.findClosestDirichletBoundary(pos, interaction))
            break;

        double radius = static_cast<double>(interaction.d) * 0.99; // slightly shrink radius to avoid numerical issues

        if (radius < epsilon) {
            auto bc = scene.boundaryCondition(interaction, pos, BoundaryType::Dirichlet);
            result += screeningThroughput * bc.a;
            return result;
        }

        // Source term: importance-sample from G(pos,·) in B(pos, dist)
        Sample<DIM> offsetSample = sampleGreensBallAtCenter<DIM>(rng, dist01, radius);
        result += screeningThroughput * scene.source(pos + offsetSample.point) *
                  greensFunctionAtCenter(offsetSample.point, radius, kappa) / offsetSample.pdf;

        // Walk to random point on sphere surface
        Sample<DIM> sphereSample = samplePoissonKernelSphereAtCenter<DIM>(rng, dist01, radius);
        pos += sphereSample.point;

        // Screening throughput update
        if (scene.isScreened()) {
            screeningThroughput *= poissonKernelAtCenter<DIM>(radius, kappa) / sphereSample.pdf;
        }
    }

    // Return NaN if maxSteps reached without hitting Dirichlet boundary
    return ScalarType::NaN();
}

// ==== Walk on Stars ====
template<typename ScalarType, int DIM>
ScalarType
walkWoSt(const PoissonScene<ScalarType, DIM>& scene,
         const Vector<DIM>& startPoint,
         int maxSteps,
         double epsilon,
         double rMin,
         pcg64& rng,
         std::uniform_real_distribution<double>& dist01,
         fcpw::Interaction<static_cast<size_t>(DIM)>& interaction)
{
    Vector<DIM> pos = startPoint;
    ScalarType result(0.0);
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;
    double screeningThroughput = 1.0;

    Vector<DIM> lastStepDir = Vector<DIM>::Zero();
    bool onNeumann = false;
    Vector<DIM> neumannNormal;
    bool flipNormal = false;

    for (int step = 0; step < maxSteps; ++step) {
        // Interior point checks
        if (scene.isInteriorOnly() && !scene.isInsideObject(pos)) {
            // Walked outside of interior domain, return NaN
            return ScalarType::NaN();
        }

        // Closest Dirichlet boundary
        if (!scene.findClosestDirichletBoundary(pos, interaction))
            break;

        double distD = static_cast<double>(interaction.d);

        if (distD < epsilon) {
            auto bc = scene.boundaryCondition(interaction, pos, BoundaryType::Dirichlet);
            result += screeningThroughput * bc.a;
            return result;
        }

        // Check if current position is nearly on Neumann boundary
        if (!onNeumann) {
            onNeumann = scene.findClosestNeumannBoundary(pos, interaction) && interaction.d < epsilon;
            flipNormal = false;
            if (onNeumann) {
                neumannNormal = interaction.n.template cast<double>().normalized();
                // Ensure normal points outward
                if (lastStepDir == Vector<DIM>::Zero()) {
                    // First step, use pos and nCheck.p to determine normal orientation
                    flipNormal = (neumannNormal.dot(interaction.p.template cast<double>() - pos) < 0.0);
                    if (flipNormal)
                        neumannNormal = -neumannNormal;
                } else {
                    flipNormal = (neumannNormal.dot(lastStepDir) < 0.0);
                    if (flipNormal)
                        neumannNormal = -neumannNormal;
                }
            }
        }

        // Closest Neumann silhouette → ball radius
        Vector<DIM> tmpPos = onNeumann ? pos - epsilon * neumannNormal
                                       : pos; // NOTE: flipNormal seems not work here, so we offset the position inward
        bool silFound = scene.findClosestNeumannSilhouette(
          tmpPos, interaction, flipNormal, 0.0f, std::numeric_limits<float>::max(), static_cast<float>(epsilon));
        double distN = silFound ? static_cast<double>(interaction.d) : std::numeric_limits<double>::max();
        double radius = std::min(distD, distN) * 0.99; // slightly shrink radius to avoid numerical issues
        radius = std::max(radius, rMin);

        // Neumann boundary contribution inside ball
        {
            if (scene.sampleNeumannBoundary(pos, static_cast<float>(radius), rng, dist01, interaction)) {
                Vector<DIM> sampleP = interaction.p.template cast<double>();
                double distSP = (sampleP - pos).norm();
                double pdf = static_cast<double>(interaction.d);

                Vector<DIM> p1 = onNeumann ? pos - epsilon * neumannNormal : pos;
                Vector<DIM> p2 = sampleP + (p1 - sampleP).normalized() * epsilon;

                if (pdf > 0.0 && distSP > EPSILON && distSP < radius && scene.checkLineOfSight(p1, p2)) {
                    double green = greensFunctionAtCenter<DIM>(sampleP - pos, radius, kappa);
                    double alpha = onNeumann ? 2.0 : 1.0;
                    auto bc = scene.boundaryCondition(interaction, p1, BoundaryType::Neumann);
                    result += screeningThroughput * alpha * bc.b * green / pdf;
                }
            }
        }

        // Sample direction
        auto [dir, pdf] = sampleUniformSphere<DIM>(rng, dist01);
        if (onNeumann) {
            if (dir.dot(neumannNormal) > 0.0)
                dir = -dir; // Flip to point inward
        }

        lastStepDir = dir;

        // Ray-intersect Neumann boundary → walk endpoint
        Vector<DIM> origin = onNeumann ? pos - epsilon * neumannNormal : pos;
        double distBoundary = radius;
        Vector<DIM> newPos;

        if (scene.rayIntersectNeumann(origin, dir, static_cast<float>(radius), interaction)) {
            onNeumann = true;
            neumannNormal = interaction.n.template cast<double>().normalized();
            flipNormal = (neumannNormal.dot(dir) < 0.0);
            if (flipNormal)
                neumannNormal = -neumannNormal; // Ensure normal points outward
            newPos = interaction.p.template cast<double>();
            distBoundary = interaction.d;
        } else {
            onNeumann = false;
            newPos = pos + dir * radius;
        }

        // Source term: importance-sample along the ray direction
        {
            auto [distSource, sourcePdf] = sampleGreensRadiusAtCenter<DIM>(rng, dist01, radius);
            if (distSource < distBoundary) {
                Vector<DIM> samplePos = pos + distSource * dir;
                ScalarType src = scene.source(samplePos);
                result += screeningThroughput * src * greensFunctionAtCenter((samplePos - pos).eval(), radius, kappa) /
                          sourcePdf;
            }
        }

        // Screening throughput update
        if (scene.isScreened()) {
            Vector<DIM> offset = newPos - pos;
            double dist = offset.norm();
            double pk = poissonKernelAtCenterOffBoundary(offset, radius, kappa);
            double areaTerm;
            if constexpr (DIM == 2) {
                areaTerm = 2.0 * PI * dist;
            } else if constexpr (DIM == 3) {
                areaTerm = 4.0 * PI * dist * dist;
            } else {
                static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
                areaTerm = 1.0;
            }
            screeningThroughput *= pk * areaTerm;
        }

        pos = newPos;
    }

    // Return NaN if maxSteps reached without hitting Dirichlet boundary
    return ScalarType::NaN();
}

// ==== Walk on Probes ====

template<typename ScalarType, int DIM>
void
walkWoP(const PoissonScene<ScalarType, DIM>& scene,
        const ProbeSet<ScalarType, DIM>& probes,
        WalkPath<ScalarType, DIM>& path,
        const Vector<DIM>& startPoint,
        int maxSteps,
        double epsilon,
        double rMin,
        pcg64& rng,
        std::uniform_real_distribution<double>& dist01,
        fcpw::Interaction<static_cast<size_t>(DIM)>& interaction)
{
    assert(maxSteps <= MaxWalkSteps);
    path.clear();

    const auto& probesVec = probes.probes();
    Vector<DIM> pos = startPoint;
    double throughput = 1.0;
    double kappa = scene.isScreened() ? scene.getScreenedKappa() : 0.0;

    // WoSt state for fallback steps
    bool onNeumann = false;
    bool flipNormal = false;
    Vector<DIM> neumannNormal;
    Vector<DIM> lastStepDir = Vector<DIM>::Zero();

    constexpr double alphaWalk = 0.7;
    constexpr double p_fb = 0.0;

    for (int step = 0; step < maxSteps; ++step) {
        // Interior point checks
        if (scene.isInteriorOnly() && !scene.isInsideObject(pos)) {
            // Walked outside of interior domain, return NaN
            return;
        }

        // ==== Dirichlet check (with probe-step pruning) ====
        bool randomFallback = (dist01(rng) < p_fb);
        bool fallbackToWoSt = onNeumann || randomFallback;

        const auto* probePtr = fallbackToWoSt ? nullptr : probes.findOptimalProbe(pos, alphaWalk);

        bool needsBoundaryCheck = (probePtr == nullptr);
        bool boundaryCheckDone = false;
        double distD = std::numeric_limits<double>::max();

        if (probePtr) {
            double gap = (1.0 - alphaWalk) * probePtr->radius;
            if (gap < epsilon)
                needsBoundaryCheck = true;
        }

        // If already decided to do WoSt step, skip boundary check. WoSt step will do its own Dirichlet check
        if (!onNeumann && needsBoundaryCheck) {
            if (scene.findClosestDirichletBoundary(pos, interaction)) {
                distD = static_cast<double>(interaction.d);
                if (distD < epsilon) {
                    auto bc = scene.boundaryCondition(interaction, pos, BoundaryType::Dirichlet);
                    path.dirichlet = throughput * bc.a;
                    path.dirichletPos = pos;
                    path.dirichletProbeIdx = (path.count > 0) ? path.probeIndices[path.count - 1] : -1;
                    return;
                }
            }
            if (scene.findClosestNeumannBoundary(pos, interaction)) {
                double distN = static_cast<double>(interaction.d);
                if (distN < epsilon) {
                    // Entered epsilon-shell of Neumann boundary, fall back to WoSt steps
                    onNeumann = true;
                    flipNormal = false;
                    if (onNeumann) {
                        neumannNormal = interaction.n.template cast<double>().normalized();
                        if (lastStepDir == Vector<DIM>::Zero()) {
                            flipNormal = (neumannNormal.dot(interaction.p.template cast<double>() - pos) < 0.0);
                            if (flipNormal)
                                neumannNormal = -neumannNormal;
                        } else {
                            flipNormal = (neumannNormal.dot(lastStepDir) < 0.0);
                            if (flipNormal)
                                neumannNormal = -neumannNormal;
                        }
                    }
                }
            }

            boundaryCheckDone = true;
        }

        fallbackToWoSt = onNeumann || (probePtr == nullptr);

        // ==== Step transition ====
        if (!fallbackToWoSt) {
            // ===== PROBE STEP =====
            double radius = probePtr->radius;
            Vector<DIM> offset = pos - probePtr->center;

            // Source term
            auto ySample = sampleGreensBall<DIM>(rng, dist01, offset, radius);
            Vector<DIM> y = probePtr->center + ySample.point;
            double green = greensFunction(offset, ySample.point, radius, kappa);
            ScalarType S_i = throughput * scene.source(y) * green / ySample.pdf;

            // Neumann: 0 (probe interior has no Neumann boundary)
            ScalarType N_i(0.0);

            // Next position z on probe boundary via Poisson kernel
            auto zSample = samplePoissonKernelSphere<DIM>(rng, dist01, offset, radius);
            Vector<DIM> newPos = probePtr->center + zSample.point;

            // Screening throughput
            if (scene.isScreened()) {
                double pk = poissonKernel<DIM>(offset, zSample.point, radius, kappa);
                throughput *= pk / zSample.pdf;
            }

            double pdfAngle = zSample.pdf;
            if constexpr (DIM == 2) {
                pdfAngle *= radius;
            } else if constexpr (DIM == 3) {
                pdfAngle *= radius * radius;
            } else {
                static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
            }
            int probeIdx = static_cast<int>(probePtr - probesVec.data());
            path.recordStep(pos, pdfAngle, S_i, N_i, probeIdx);
            lastStepDir = (newPos - pos).normalized();
            pos = newPos;
            onNeumann = false;
        } else {
            // ===== FALLBACK WoSt STEP =====
            // Ball radius from closest boundaries

            // Dirichlet check (if not done by probe pruning)
            if (!boundaryCheckDone) {
                if (scene.findClosestDirichletBoundary(pos, interaction)) {
                    distD = static_cast<double>(interaction.d);
                }
            }

            double distN = std::numeric_limits<double>::max();
            Vector<DIM> tmpPos = onNeumann ? pos - epsilon * neumannNormal : pos;
            if (scene.findClosestNeumannSilhouette(tmpPos,
                                                   interaction,
                                                   flipNormal,
                                                   0.0f,
                                                   std::numeric_limits<float>::max(),
                                                   static_cast<float>(epsilon)))
                distN = static_cast<double>(interaction.d);

            double radius = std::min(distD, distN) * 0.99; // slightly shrink radius to avoid numerical issues
            radius = std::max(radius, rMin);

            // Neumann boundary contribution inside ball
            ScalarType N_i(0.0);
            {
                if (scene.sampleNeumannBoundary(pos, static_cast<float>(radius), rng, dist01, interaction)) {
                    Vector<DIM> sampleP = interaction.p.template cast<double>();
                    double ndistSP = (sampleP - pos).norm();
                    double nPdf = static_cast<double>(interaction.d);

                    Vector<DIM> p1 = onNeumann ? pos - epsilon * neumannNormal : pos;
                    Vector<DIM> p2 = sampleP + (p1 - sampleP).normalized() * epsilon;

                    if (nPdf > 0.0 && ndistSP > EPSILON && ndistSP < radius && scene.checkLineOfSight(p1, p2)) {
                        double green = greensFunctionAtCenter<DIM>(sampleP - pos, radius, kappa);
                        double alpha = onNeumann ? 2.0 : 1.0;
                        auto bc = scene.boundaryCondition(interaction, p1, BoundaryType::Neumann);
                        N_i = throughput * alpha * bc.b * green / nPdf;
                    }
                }
            }

            // Sample direction
            auto [dir, dirPdf] = sampleUniformSphere<DIM>(rng, dist01);
            if (onNeumann) {
                if (dir.dot(neumannNormal) > 0.0)
                    dir = -dir;
            }
            lastStepDir = dir;

            // Ray-intersect Neumann → walk endpoint
            Vector<DIM> origin = onNeumann ? pos - epsilon * neumannNormal : pos;
            double distBoundary = radius;
            Vector<DIM> newPos;

            if (scene.rayIntersectNeumann(origin, dir, static_cast<float>(radius), interaction)) {
                onNeumann = true;
                neumannNormal = interaction.n.template cast<double>().normalized();
                flipNormal = (neumannNormal.dot(dir) < 0.0);
                if (flipNormal)
                    neumannNormal = -neumannNormal;
                newPos = interaction.p.template cast<double>();
                distBoundary = interaction.d;
            } else {
                onNeumann = false;
                newPos = pos + dir * radius;
                distBoundary = radius;
            }

            // Source term along ray direction
            ScalarType S_i(0.0);
            {
                auto [distSource, sourcePdf] = sampleGreensRadiusAtCenter<DIM>(rng, dist01, radius);
                if (distSource < distBoundary) {
                    Vector<DIM> samplePos = pos + distSource * dir;
                    S_i = throughput * scene.source(samplePos) *
                          greensFunctionAtCenter<DIM>((samplePos - pos).eval(), radius, kappa) / sourcePdf;
                }
            }

            // Screening throughput
            if (scene.isScreened()) {
                Vector<DIM> offset = newPos - pos;
                double dist = offset.norm();
                double pk = poissonKernelAtCenterOffBoundary<DIM>(offset, radius, kappa);
                double areaTerm;
                if constexpr (DIM == 2)
                    areaTerm = 2.0 * PI * dist;
                else if constexpr (DIM == 3)
                    areaTerm = 4.0 * PI * dist * dist;
                else
                    areaTerm = 1.0;
                throughput *= pk * areaTerm;
            }

            path.recordStep(pos, 0.5 * INV_PI, S_i, N_i, -1);
            pos = newPos;
        }
    }
    return;
}

WOS_NAMESPACE_CLOSE_SCOPE
