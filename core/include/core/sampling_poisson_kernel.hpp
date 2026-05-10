#pragma once

#include "core/math_defs.hpp"
#include "core/sample.hpp"
#include "core/sampling_uniform.hpp"
#include "core/vector.hpp"

#include "pcg_random.hpp"

#include <random>

WOS_NAMESPACE_OPEN_SCOPE

// TODO: Importance sample y ∈ ∂B(0,R) according to Poisson kernel P(x, ·)
// Returns sample point and its pdf

template<int DIM>
[[nodiscard]] Sample<DIM>
samplePoissonKernelSphereAtCenter(pcg64& rng, std::uniform_real_distribution<double>& dist01, double R)
{
    Sample<DIM> s = sampleUniformSphere<DIM>(rng, dist01);
    s.point *= R;
    if constexpr (DIM == 2) {
        s.pdf /= R;
    } else if constexpr (DIM == 3) {
        s.pdf /= (R * R);
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
    return s;
}

template<int DIM>
[[nodiscard]] Sample<DIM>
samplePoissonKernelSphereOffCenter(pcg64& rng,
                                   std::uniform_real_distribution<double>& dist01,
                                   const Vector<DIM>& x,
                                   double R)
{
    if constexpr (DIM == 2) {
        double rho = x.norm();

        // Wrapped Cauchy inverse transform
        double phi0 = std::atan2(x.y(), x.x());
        double factor = (R - rho) / (R + rho);
        double u = dist01(rng);
        double theta_rel = 2.0 * std::atan(factor * std::tan(PI * (u - 0.5)));
        double theta = theta_rel + phi0;
        Vector<DIM> point{ R * std::cos(theta), R * std::sin(theta) };

        double d_sq = (x - point).squaredNorm();
        double pdf = 0.5 * INV_PI * (R * R - rho * rho) / (R * d_sq);
        return { point, pdf };
    } else if constexpr (DIM == 3) {
        double rho = x.norm();

        // Inverse transform: 1/D uniform in [1/(R+rho), 1/(R-rho)]
        double v_min = 1.0 / (R + rho);
        double v_max = 1.0 / (R - rho);
        double u1 = dist01(rng);
        double inv_D = v_min + u1 * (v_max - v_min);
        double D = 1.0 / std::max(inv_D, EPSILON);

        // cos(theta) via law of cosines
        double R_sq = R * R;
        double rho_sq = rho * rho;
        double mu = (R_sq + rho_sq - D * D) / (2.0 * R * rho);
        double mu_clamped = std::clamp(mu, -1.0, 1.0);
        double sinT = std::sqrt(std::max(0.0, 1.0 - mu_clamped * mu_clamped));
        double phi = 2.0 * PI * dist01(rng);

        // Build point in local frame (z-axis aligned with x)
        Vector<DIM> local{ R * sinT * std::cos(phi), R * sinT * std::sin(phi), R * mu_clamped };

        // Rotate (0,0,1) to x/rho via Rodrigues formula
        Vector<DIM> target = x / rho;
        Vector<DIM> zAxis{ 0.0, 0.0, 1.0 };
        Vector<DIM> v = zAxis.cross(target);
        double c = zAxis.dot(target);
        double s_sq = v.squaredNorm();

        Vector<DIM> point;
        if (s_sq < EPSILON_SQ) {
            point = (c > 0.0) ? local : Vector<DIM>{ -local.x(), -local.y(), -local.z() };
        } else {
            // Rodrigues: R = I + K + K²·(1-c)/s²  →  R·v = v + K×v + ((1-c)/s²)·K×(K×v)
            Vector<DIM> Kv = v.cross(local);
            Vector<DIM> KKv = ((1.0 - c) / s_sq) * v.cross(Kv);
            point = local + Kv + KKv;
        }

        double d_sq = (x - point).squaredNorm();
        double d = std::sqrt(d_sq);
        double pdf = 0.25 * INV_PI * (R_sq - rho_sq) / (R * d_sq * d);
        return { point, pdf };
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
    }
}

template<int DIM>
[[nodiscard]] Sample<DIM>
samplePoissonKernelSphere(pcg64& rng, std::uniform_real_distribution<double>& dist01, const Vector<DIM>& x, double R)
{
    double x_squared_norm = x.squaredNorm();
    if (x_squared_norm < EPSILON_SQ) {
        return samplePoissonKernelSphereAtCenter<DIM>(rng, dist01, R);
    } else {
        return samplePoissonKernelSphereOffCenter<DIM>(rng, dist01, x, R);
    }
}

WOS_NAMESPACE_CLOSE_SCOPE
