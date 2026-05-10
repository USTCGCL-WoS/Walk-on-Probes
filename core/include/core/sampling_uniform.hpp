#pragma once

#include "core/math_defs.hpp"
#include "core/sample.hpp"
#include "core/vector.hpp"

#include "pcg_random.hpp"

#include <cmath>
#include <random>

WOS_NAMESPACE_OPEN_SCOPE

template<int DIM>
[[nodiscard]] Sample<DIM>
sampleUniformSphere(pcg64& rng, std::uniform_real_distribution<double>& dist01)
{
    if constexpr (DIM == 2) {
        double theta = dist01(rng) * 2.0 * PI;
        return { { std::cos(theta), std::sin(theta) }, 0.5 * INV_PI };
    } else {
        double z = 2.0 * dist01(rng) - 1.0;
        double theta = dist01(rng) * 2.0 * PI;
        double r = std::sqrt(1.0 - z * z);
        return { { r * std::cos(theta), r * std::sin(theta), z }, 0.25 * INV_PI };
    }
}

template<int DIM>
[[nodiscard]] Sample<DIM>
sampleUniformBall(pcg64& rng, std::uniform_real_distribution<double>& dist01, double R)
{
    if constexpr (DIM == 2) {
        double r = R * std::sqrt(dist01(rng));
        double theta = dist01(rng) * 2.0 * PI;
        double pdf = INV_PI / (R * R);
        return { { r * std::cos(theta), r * std::sin(theta) }, pdf };
    } else {
        double r = R * std::cbrt(dist01(rng));
        double z = 2.0 * dist01(rng) - 1.0;
        double theta = dist01(rng) * 2.0 * PI;
        double rxy = r * std::sqrt(1.0 - z * z);
        double pdf = 3.0 / (4.0 * PI * R * R * R);
        return { { rxy * std::cos(theta), rxy * std::sin(theta), r * z }, pdf };
    }
}

WOS_NAMESPACE_CLOSE_SCOPE
