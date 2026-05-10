#pragma once

#include "core/math_defs.hpp"
#include "core/vector.hpp"
#include <cassert>
#include <cmath>

WOS_NAMESPACE_OPEN_SCOPE

// TODO: Green's function for screened Poisson equation (Δ - κ²)G = δ in ball B(0,R)

template<int DIM>
[[nodiscard]] double
greensFunctionAtCenter(const Vector<DIM>& y, double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        double r = y.norm();
        if (std::abs(kappa) < EPSILON) {
            return 0.5 * INV_PI * std::log(R / r);
        } else {
            double kr = kappa * r;
            double kR = kappa * R;
            return 0.5 * INV_PI *
                   (std::cyl_bessel_k(0, kr) -
                    std::cyl_bessel_k(0, kR) * std::cyl_bessel_i(0, kr) / std::cyl_bessel_i(0, kR));
        }
    } else if constexpr (DIM == 3) {
        double r = y.norm();
        if (std::abs(kappa) < EPSILON) {
            return 0.25 * INV_PI * (1.0 / r - 1.0 / R);
        } else {
            double kr = kappa * r;
            double kR = kappa * R;
            return 0.25 * INV_PI * ((std::exp(-kr) / r) - (std::exp(-kR) * std::sinh(kr) / (r * std::sinh(kR))));
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0;
    }
}

template<int DIM>
[[nodiscard]] double
greensFunctionOffCenter(const Vector<DIM>& x, const Vector<DIM>& y, double R, double kappa = 0.0) noexcept
{
    // TODO
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            Vector<DIM> x_asterisk = x * R * R / x.squaredNorm();
            double r = (y - x).norm();
            double r_asterisk = (y - x_asterisk).norm();
            return 0.5 * INV_PI * std::log(x.norm() * r_asterisk / (R * r));
        } else {
            double x_norm = x.norm();
            double y_norm = y.norm();
            double r_minus = std::min(x_norm, y_norm);
            double r_plus = std::max(x_norm, y_norm);
            double theta = std::acos(std::clamp(x.dot(y) / (x_norm * y_norm), -1.0, 1.0));
            double kr_minus = kappa * r_minus;
            double kr_plus = kappa * r_plus;
            double kR = kappa * R;

            double sum = 0.0;
            double bessel_term = 0.0;
            double cos_term = 0.0;
            int n = 0;
            constexpr int max_iter = 1000;
            for (n = 0; n < max_iter; ++n) {
                bessel_term = (n == 0 ? 1.0 : 2.0) * std::cyl_bessel_i(n, kr_minus) *
                              (std::cyl_bessel_k(n, kr_plus) -
                               std::cyl_bessel_k(n, kR) * std::cyl_bessel_i(n, kr_plus) / std::cyl_bessel_i(n, kR));
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
            Vector<DIM> x_asterisk = x * R * R / x.squaredNorm();
            double r = (y - x).norm();
            double r_asterisk = (y - x_asterisk).norm();
            return 0.25 * INV_PI * (1.0 / r - R / (x.norm() * r_asterisk));
        } else {
            double x_norm = x.norm();
            double y_norm = y.norm();
            double r_minus = std::min(x_norm, y_norm);
            double r_plus = std::max(x_norm, y_norm);
            double cos_theta = x.dot(y) / (x_norm * y_norm);
            double kr_minus = kappa * r_minus;
            double kr_plus = kappa * r_plus;
            double kR = kappa * R;

            double sum = 0.0;
            double bessel_term = 0.0;
            // Legendre recurrence: P_0=1, P_1=x, P_{n}=((2n-1)*x*P_{n-1}-(n-1)*P_{n-2})/n
            double P_nm2 = 1.0;
            double P_nm1 = 1.0;
            int n = 0;
            constexpr int max_iter = 1000;
            for (n = 0; n < max_iter; ++n) {
                double P_n;
                if (n == 0) {
                    P_n = 1.0;
                } else if (n == 1) {
                    P_n = cos_theta;
                } else {
                    P_n = ((2.0 * n - 1.0) * cos_theta * P_nm1 - (n - 1.0) * P_nm2) / n;
                }

                bessel_term = (2.0 * n + 1.0) * std::cyl_bessel_i(n + 0.5, kr_minus) *
                              (std::cyl_bessel_k(n + 0.5, kr_plus) - std::cyl_bessel_k(n + 0.5, kR) *
                                                                       std::cyl_bessel_i(n + 0.5, kr_plus) /
                                                                       std::cyl_bessel_i(n + 0.5, kR)) /
                              std::sqrt(r_minus * r_plus);
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
greensFunction(const Vector<DIM>& x, const Vector<DIM>& y, double R, double kappa = 0.0) noexcept
{
    double x_squared_norm = x.squaredNorm();
    double y_squared_norm = y.squaredNorm();
    if (x_squared_norm < EPSILON_SQ) {
        return greensFunctionAtCenter(y, R, kappa);
    } else if (y_squared_norm < EPSILON_SQ) {
        return greensFunctionAtCenter(x, R, kappa);
    } else {
        return greensFunctionOffCenter(x, y, R, kappa);
    }
}

template<int DIM>
[[nodiscard]] double
greensNormalizationAtCenter(double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            return 0.25 * R * R;
        } else {
            double kR = kappa * R;
            return (1.0 - 1.0 / std::cyl_bessel_i(0, kR)) / (kappa * kappa);
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            return R * R / 6.0;
        } else {
            double kR = kappa * R;
            return (1.0 - (kR / std::sinh(kR))) / (kappa * kappa);
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0;
    }
}

template<int DIM>
[[nodiscard]] double
greensNormalizationOffCenter(const Vector<DIM>& x, double R, double kappa = 0.0) noexcept
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            return 0.25 * (R * R - x.squaredNorm());
        } else {
            // TODO: screened Poisson off-center normalization
            assert(false && "Screened Poisson off-center normalization not implemented");
            return 0.0;
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            return (R * R - x.squaredNorm()) / 6.0;
        } else {
            // TODO: screened Poisson off-center normalization
            assert(false && "Screened Poisson off-center normalization not implemented");
            return 0.0;
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0;
    }
}

template<int DIM>
[[nodiscard]] double
greensNormalization(const Vector<DIM>& x, double R, double kappa = 0.0) noexcept
{
    double x_squared_norm = x.squaredNorm();
    if (x_squared_norm < EPSILON_SQ) {
        return greensNormalizationAtCenter<DIM>(R, kappa);
    } else {
        return greensNormalizationOffCenter<DIM>(x, R, kappa);
    }
}

WOS_NAMESPACE_CLOSE_SCOPE
