#pragma once

#include "core/green.hpp"
#include "core/math_defs.hpp"
#include "core/sample.hpp"
#include "core/vector.hpp"

#include "pcg_random.hpp"

#include <random>

WOS_NAMESPACE_OPEN_SCOPE

// TODO: Importance sample y ∈ B(0,R) according to Green's function G(x, ·)
// Returns sample point and its pdf

template<int DIM>
[[nodiscard]] SampleLength<DIM>
sampleGreensRadiusAtCenter(pcg64& rng, std::uniform_real_distribution<double>& dist01, double R)
{
    if constexpr (DIM == 2) {
        double u1 = 1.0 - dist01(rng);
        double u2 = 1.0 - dist01(rng);
        double r = R * std::sqrt(u1 * u2);
        double pdf = 2.0 * INV_PI * std::log(R / r) / (R * R); // Directly return PDF of area
        return { r, pdf };
    } else if constexpr (DIM == 3) {
        // r/R ~ Beta(2,2) via median-of-3 uniforms
        double u1 = 1.0 - dist01(rng);
        double u2 = 1.0 - dist01(rng);
        double u3 = 1.0 - dist01(rng);
        if (u1 > u2)
            std::swap(u1, u2);
        if (u2 > u3)
            std::swap(u2, u3);
        if (u1 > u2)
            std::swap(u1, u2);
        double t = u2;
        double r = R * t;
        double pdf = 1.5 * INV_PI * (1.0 / r - 1.0 / R) / (R * R); // Directly return PDF of volume
        return { r, pdf };
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
}

template<int DIM>
[[nodiscard]] Sample<DIM>
sampleGreensBallAtCenter(pcg64& rng, std::uniform_real_distribution<double>& dist01, double R)
{
    if constexpr (DIM == 2) {
        auto [r, pdf] = sampleGreensRadiusAtCenter<DIM>(rng, dist01, R);
        double theta = 2.0 * PI * dist01(rng);
        Vector<DIM> point{ r * std::cos(theta), r * std::sin(theta) };
        return { point, pdf };
    } else if constexpr (DIM == 3) {
        auto [r, pdf] = sampleGreensRadiusAtCenter<DIM>(rng, dist01, R);
        double phi = std::acos(1.0 - 2.0 * dist01(rng)); // cos(phi) ~ Uniform(-1,1)
        double theta = 2.0 * PI * dist01(rng);
        Vector<DIM> point{ r * std::sin(phi) * std::cos(theta),
                           r * std::sin(phi) * std::sin(theta),
                           r * std::cos(phi) };
        return { point, pdf };
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
}

template<int DIM>
[[nodiscard]] Sample<DIM>
sampleGreensBallOffCenter(pcg64& rng, std::uniform_real_distribution<double>& dist01, const Vector<DIM>& x, double R)
{
    // Rejection sampling via bounding sphere B(x, R+|x|) ⊇ B(0,R)
    // Proposal: at-center Green's function G(0,·) in the bounding sphere

    if constexpr (DIM == 2 || DIM == 3) {
        double R_b = R + x.norm();
        double R_sq = R * R;
        double G_target = 0.0;
        double G_proposal = 0.0;
        double G_target_norm = greensNormalizationOffCenter(x, R);
        double P_accept = 0.0;
        constexpr int max_iter = 1000;
        for (int iter = 0; iter < max_iter; ++iter) {
            auto [y_local, pdf] = sampleGreensBallAtCenter<DIM>(rng, dist01, R_b);
            Vector<DIM> y = x + y_local;
            if (y.squaredNorm() > R_sq) {
                continue;
            }
            G_target = greensFunctionOffCenter(x, y, R);
            G_proposal = greensFunctionAtCenter(y_local, R_b);
            P_accept = G_target / G_proposal;
            if (dist01(rng) < P_accept) {
                return { y, G_target / G_target_norm };
            }
        }
        return { x, 0.0 };
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
}

template<int DIM>
[[nodiscard]] Sample<DIM>
sampleGreensBall(pcg64& rng, std::uniform_real_distribution<double>& dist01, const Vector<DIM>& x, double R)
{
    double x_squared_norm = x.squaredNorm();
    if (x_squared_norm < EPSILON_SQ) {
        return sampleGreensBallAtCenter<DIM>(rng, dist01, R);
    } else {
        return sampleGreensBallOffCenter<DIM>(rng, dist01, x, R);
    }
}

WOS_NAMESPACE_CLOSE_SCOPE
