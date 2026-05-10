#pragma once

#include "core/math_defs.hpp"
#include "core/vector.hpp"

WOS_NAMESPACE_OPEN_SCOPE

// TODO: Poisson kernel for screened Poisson equation (Δ - κ²)u = f in ball B(0,R)

template<int DIM>
[[nodiscard]] double
poissonKernelAtCenter(double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            return 0.5 * INV_PI / R;
        } else {
            double kR = kappa * R;
            return 0.5 * INV_PI / (R * std::cyl_bessel_i(0, kR));
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            return 0.25 * INV_PI / (R * R);
        } else {
            double kR = kappa * R;
            return 0.25 * INV_PI * kR / (R * R * std::sinh(kR));
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0; // Unreachable
    }
}

template<int DIM>
[[nodiscard]] double
poissonKernelAtCenterOffBoundary(const Vector<DIM>& z, double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            double r = z.norm();
            return 0.5 * INV_PI / r;
        } else {
            double r = z.norm();
            double kR = kappa * R;
            double kr = kappa * r;
            return 0.5 * INV_PI * kappa *
                   (std::cyl_bessel_k(1, kr) +
                    std::cyl_bessel_i(1, kr) * std::cyl_bessel_k(0, kR) / std::cyl_bessel_i(0, kR));
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            double r = z.norm();
            return 0.25 * INV_PI / (r * r);
        } else {
            double r = z.norm();
            double kR = kappa * R;
            double kr = kappa * r;
            return 0.25 * INV_PI *
                   (std::exp(-kr) * (kr + 1) + std::exp(-kR) * (std::cosh(kr) * kr - std::sinh(kr)) / std::sinh(kR)) /
                   (r * r);
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0; // Unreachable
    }
}

template<int DIM>
[[nodiscard]] double
poissonKernelOffCenter(const Vector<DIM>& x, const Vector<DIM>& z, double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            double r = x.norm();
            double d = (x - z).norm();
            return 0.5 * INV_PI * (R * R - r * r) / (R * d * d);
        } else {
            double r = x.norm();
            double theta = std::acos(std::clamp(x.dot(z) / (r * R), -1.0, 1.0));
            double kR = kappa * R;
            double kr = kappa * r;

            double sum = 0.0;
            double bessel_term = 0.0;
            double cos_term = 0.0;
            int n = 0;
            constexpr int max_iter = 200;
            for (n = 0; n < max_iter; ++n) {
                bessel_term = (n == 0 ? 1.0 : 2.0) * std::cyl_bessel_i(n, kr) / (std::cyl_bessel_i(n, kR) * R);
                cos_term = std::cos(n * theta);
                sum += bessel_term * cos_term;
                if (bessel_term < EPSILON) {
                    break;
                }
            }
            return 0.5 * INV_PI * sum;
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            double r = x.norm();
            double d = (x - z).norm();
            return 0.25 * INV_PI * (R * R - r * r) / (R * d * d * d);
        } else {
            double r = x.norm();
            double cos_theta = x.dot(z) / (r * R);
            double kR = kappa * R;
            double kr = kappa * r;

            double sum = 0.0;
            double bessel_term = 0.0;
            // Legendre recurrence: P_0=1, P_1=x, P_{n}=((2n-1)*x*P_{n-1}-(n-1)*P_{n-2})/n
            double P_nm2 = 1.0;
            double P_nm1 = 1.0;
            int n = 0;
            constexpr int max_iter = 200;
            for (n = 0; n < max_iter; ++n) {
                double P_n;
                if (n == 0) {
                    P_n = 1.0;
                } else if (n == 1) {
                    P_n = cos_theta;
                } else {
                    P_n = ((2.0 * n - 1.0) * cos_theta * P_nm1 - (n - 1.0) * P_nm2) / n;
                }

                bessel_term = (2.0 * n + 1.0) * std::cyl_bessel_i(n + 0.5, kr) /
                              (std::cyl_bessel_i(n + 0.5, kR) * R * std::sqrt(R * r));
                sum += bessel_term * P_n;
                if (bessel_term < EPSILON) {
                    break;
                }

                P_nm2 = P_nm1;
                P_nm1 = P_n;
            }
            return 0.25 * INV_PI * sum;
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0;
    }
}

template<int DIM>
[[nodiscard]] double
poissonKernel(const Vector<DIM>& x, const Vector<DIM>& z, double R, double kappa = 0.0) noexcept
{
    double x_squared_norm = x.squaredNorm();
    if (x_squared_norm < EPSILON_SQ) {
        return poissonKernelAtCenter<DIM>(R, kappa);
    } else {
        return poissonKernelOffCenter<DIM>(x, z, R, kappa);
    }
}

WOS_NAMESPACE_CLOSE_SCOPE
