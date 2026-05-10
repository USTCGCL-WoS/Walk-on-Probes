#pragma once

#include <cmath>
#include <mutex>
#include <vector>

#include "core/math_defs.hpp"
#include "core/scalar.hpp"
#include "core/vector.hpp"

#include <spdlog/spdlog.h>

WOS_NAMESPACE_OPEN_SCOPE

// Movable/copyable mutex wrapper — the wrapped std::mutex is non-movable/non-copyable,
// so copy/move of this wrapper creates a fresh mutex. This lets Probe stay copyable/movable
// while each copy gets its own independent mutex.
struct MovableMutex
{
    std::mutex mtx;
    MovableMutex() = default;
    MovableMutex(MovableMutex&&) noexcept {}
    MovableMutex& operator=(MovableMutex&&) noexcept { return *this; }
    MovableMutex(const MovableMutex&) noexcept {}
    MovableMutex& operator=(const MovableMutex&) noexcept { return *this; }
    void lock() { mtx.lock(); }
    void unlock() { mtx.unlock(); }
};

// ==== radialBasis ====

template<int DIM>
double
radialBasis(double r, double R, int l, double kappa)
{
    if constexpr (DIM == 2) {
        if (std::abs(kappa) < EPSILON) {
            if (r < EPSILON)
                return (l == 0) ? 1.0 : 0.0;
            return std::pow(r / R, static_cast<double>(l));
        } else {
            double kr = kappa * r;
            double kR = kappa * R;
            return std::cyl_bessel_i(l, kr) / std::cyl_bessel_i(l, kR);
        }
    } else if constexpr (DIM == 3) {
        if (std::abs(kappa) < EPSILON) {
            if (r < EPSILON)
                return (l == 0) ? 1.0 : 0.0;
            return std::pow(r / R, static_cast<double>(l));
        } else {
            if (r < EPSILON)
                return (l == 0) ? kappa * R / std::sinh(kappa * R) : 0.0;
            double kr = kappa * r;
            double kR = kappa * R;
            return std::cyl_bessel_i(l + 0.5, kr) / std::cyl_bessel_i(l + 0.5, kR) * std::sqrt(kR / kr);
        }
    } else {
        static_assert(DIM == 2 || DIM == 3, "Unsupported dimension");
        return 0.0;
    }
}

// Forward declaration
template<typename ScalarType, int DIM>
struct Probe;

// ==== DIM=2 specialization — Fourier series expansion ====

template<typename ScalarType>
struct Probe<ScalarType, 2>
{
    Vector<2> center;
    double radius = 0.0;

    // Fourier accumulators — per-sample online accumulation
    std::vector<ScalarType> sumWUcos; // l=0..L, Σ w_i * u_i * cos(lθ_i)
    std::vector<ScalarType> sumWUsin; // l=0..L, Σ w_i * u_i * sin(lθ_i)
    std::vector<double> sumWcos;      // l=0..L, Σ w_i * cos(lθ_i)
    std::vector<double> sumWsin;      // l=0..L, Σ w_i * sin(lθ_i)
    double weightSum = 0.0;           // W_Σ = Σ w_i
    int sampleCount = 0;              // total raw sample count

    // Finalized coefficients
    int L = 0;                       // number of Fourier modes
    std::vector<ScalarType> aCoeffs; // l=0..L
    std::vector<ScalarType> bCoeffs; // l=0..L (b[0]=0)

    // Source accumulation per cached target point
    std::vector<ScalarType> srcAccum;

    // Cached target point indices (into solver targetPoints_)
    std::vector<int> targetIndices;

    int srcCount = 0;

    MovableMutex mutex_;

    bool contains(const Vector<2>& p, double alpha) const
    {
        return (p - center).squaredNorm() <= radius * radius * alpha * alpha;
    }

    void initExpansion(int L)
    {
        this->L = L;
        sumWUcos.assign(L + 1, ScalarType(0.0));
        sumWUsin.assign(L + 1, ScalarType(0.0));
        sumWcos.assign(L + 1, 0.0);
        sumWsin.assign(L + 1, 0.0);
        weightSum = 0.0;
        sampleCount = 0;
        aCoeffs.assign(L + 1, ScalarType(0.0));
        bCoeffs.assign(L + 1, ScalarType(0.0));
    }

    void addSample(Vector<2> samplePos, const ScalarType& u, double pdf)
    {
        std::lock_guard lock(mutex_);
        addSampleUnlocked(samplePos, u, pdf);
    }

    void addSampleUnlocked(Vector<2> samplePos, const ScalarType& u, double pdf)
    {
        double w = 0.5 * INV_PI / pdf;
        double angle = std::atan2(samplePos[1] - center[1], samplePos[0] - center[0]);
        double cosLtheta, sinLtheta;

        ScalarType wu = w * u;
        for (int l = 0; l <= L; ++l) {
            cosLtheta = std::cos(l * angle);
            sinLtheta = std::sin(l * angle);
            sumWUcos[l] += wu * cosLtheta;
            sumWUsin[l] += wu * sinLtheta;
            sumWcos[l] += w * cosLtheta;
            sumWsin[l] += w * sinLtheta;
        }
        weightSum += w;
        sampleCount += 1;
    }

    void finalizeCoefficients(bool useSN, bool useCV)
    {
        if (sampleCount == 0)
            return;
        double denom = useSN ? weightSum : static_cast<double>(sampleCount);

        ScalarType a0 = sumWUcos[0] / denom; // cos(0*θ)=1
        aCoeffs[0] = a0;
        bCoeffs[0] = ScalarType(0.0);

        for (size_t l = 1; l < aCoeffs.size(); ++l) {
            ScalarType mean_u_cos = sumWUcos[l] / denom;
            ScalarType mean_u_sin = sumWUsin[l] / denom;
            if (useCV) {
                double mean_cos = sumWcos[l] / denom;
                double mean_sin = sumWsin[l] / denom;
                aCoeffs[l] = mean_u_cos - a0 * mean_cos;
                bCoeffs[l] = mean_u_sin - a0 * mean_sin;
            } else {
                aCoeffs[l] = mean_u_cos;
                bCoeffs[l] = mean_u_sin;
            }
        }

        if (srcCount > 0)
            for (auto& s : srcAccum)
                s /= srcCount;
    }

    ScalarType evaluate(const Vector<2>& p, double kappa) const
    {
        Vector<2> offset = p - center;
        double d = offset.norm();
        if (d < EPSILON)
            return aCoeffs[0] * radialBasis<2>(0.0, radius, 0, kappa);
        else {
            double theta = std::atan2(offset[1], offset[0]);

            ScalarType val = aCoeffs[0] * radialBasis<2>(d, radius, 0, kappa);
            for (size_t l = 1; l < aCoeffs.size(); ++l) {
                double Rl = radialBasis<2>(d, radius, static_cast<int>(l), kappa);
                val += 2.0 * Rl * (aCoeffs[l] * std::cos(l * theta) + bCoeffs[l] * std::sin(l * theta));
            }
            if (val.isNaN()) {
                return ScalarType::NaN();
            }
            return val;
        }
    }
};

// ==== DIM=3 specialization — spherical harmonic expansion ====

template<typename ScalarType>
struct Probe<ScalarType, 3>
{
    Vector<3> center;
    double radius = 0.0;

    // Triangular storage: shIdx(l,m) = l*(l+1)/2 + m,  m ∈ [0, l],  total = (L+1)*(L+2)/2
    int L = 0;
    int nCoeff = 0;

    // Per-sample spherical harmonic accumulators
    // Y_l^m ∝ P_l^m(cosθ)·cos(mφ) stored in cosM (m ≥ 0)
    // Y_l^{-m} ∝ P_l^m(cosθ)·sin(mφ) stored in sinM (m > 0)
    std::vector<ScalarType> sumWU_cosM; // Σ w*u*Y_l^m  (m ≥ 0, cos branch)
    std::vector<ScalarType> sumWU_sinM; // Σ w*u*Y_l^{-m} (m > 0, sin branch)
    std::vector<double> sumW_cosM;      // Σ w*Y_l^m  (m ≥ 0, for CV)
    std::vector<double> sumW_sinM;      // Σ w*Y_l^{-m} (m > 0, for CV)

    // Finalized coefficients
    std::vector<ScalarType> aCoeffs; // a_l^m (cos branch, m ≥ 0)
    std::vector<ScalarType> bCoeffs; // b_l^m (sin branch, m > 0, b_l^0 unused)

    // Precomputed normalization constants N_l^m (same triangular layout)
    std::vector<double> normFactors;

    double weightSum = 0.0;
    int sampleCount = 0;

    // Source accumulation per cached target point
    std::vector<ScalarType> srcAccum;
    std::vector<int> targetIndices;
    int srcCount = 0;

    MovableMutex mutex_;

    static constexpr int shIdx(int l, int m) { return l * (l + 1) / 2 + m; }

    bool contains(const Vector<3>& p, double alpha) const
    {
        return (p - center).squaredNorm() <= radius * radius * alpha * alpha;
    }

    void initExpansion(int L)
    {
        this->L = L;
        nCoeff = (L + 1) * (L + 2) / 2;
        sumWU_cosM.assign(nCoeff, ScalarType(0.0));
        sumWU_sinM.assign(nCoeff, ScalarType(0.0));
        sumW_cosM.assign(nCoeff, 0.0);
        sumW_sinM.assign(nCoeff, 0.0);
        weightSum = 0.0;
        sampleCount = 0;
        aCoeffs.assign(nCoeff, ScalarType(0.0));
        bCoeffs.assign(nCoeff, ScalarType(0.0));

        // Precompute normalization constants N_l^m
        normFactors.assign(nCoeff, 0.0);
        for (int l = 0; l <= L; ++l) {
            for (int m = 0; m <= l; ++m) {
                int idx = shIdx(l, m);
                // K_l^m = sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
                // Using lgamma for stable computation of factorial ratio
                double lnK = 0.5 * (std::log(2.0 * l + 1.0) - std::log(4.0 * PI) + std::lgamma(l - m + 1.0) -
                                    std::lgamma(l + m + 1.0));
                double K = std::exp(lnK);
                if (m > 0)
                    K *= SQRT2;
                normFactors[idx] = K;
            }
        }
    }

    void addSample(Vector<3> samplePos, const ScalarType& u, double pdf)
    {
        std::lock_guard lock(mutex_);
        addSampleUnlocked(samplePos, u, pdf);
    }

    void addSampleUnlocked(Vector<3> samplePos, const ScalarType& u, double pdf)
    {
        double w = 0.25 * INV_PI / pdf;

        Vector<3> offset = samplePos - center;
        double r = offset.norm();
        if (r < EPSILON)
            return;

        double cosTheta = offset[2] / r;
        double phi = std::atan2(offset[1], offset[0]);

        // Precompute cos(mφ) and sin(mφ) via Chebyshev recurrence
        std::vector<double> cosMPhi(L + 1), sinMPhi(L + 1);
        cosMPhi[0] = 1.0;
        sinMPhi[0] = 0.0;
        if (L >= 1) {
            double cosPhi = std::cos(phi);
            double sinPhi = std::sin(phi);
            cosMPhi[1] = cosPhi;
            sinMPhi[1] = sinPhi;
            for (int m = 2; m <= L; ++m) {
                cosMPhi[m] = 2.0 * cosPhi * cosMPhi[m - 1] - cosMPhi[m - 2];
                sinMPhi[m] = 2.0 * cosPhi * sinMPhi[m - 1] - sinMPhi[m - 2];
            }
        }

        // Precompute associated Legendre P_l^m(cosθ) for all l=0..L, m=0..l
        // Without Condon-Shortley phase
        std::vector<double> P(nCoeff, 0.0);
        double sx = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta)); // sinθ
        P[shIdx(0, 0)] = 1.0;

        // Diagonal: P_m^m
        for (int m = 1; m <= L; ++m) {
            P[shIdx(m, m)] = (2.0 * m - 1.0) * sx * P[shIdx(m - 1, m - 1)];
        }

        // Off-diagonal: P_l^m for l > m
        for (int m = 0; m <= L; ++m) {
            if (m + 1 <= L) {
                P[shIdx(m + 1, m)] = (2.0 * m + 1.0) * cosTheta * P[shIdx(m, m)];
            }
            for (int l = m + 2; l <= L; ++l) {
                double term1 = (2.0 * l - 1.0) * cosTheta * P[shIdx(l - 1, m)];
                double term2 = (l + m - 1.0) * P[shIdx(l - 2, m)];
                P[shIdx(l, m)] = (term1 - term2) / (l - m);
            }
        }

        ScalarType wu = w * u;
        for (int l = 0; l <= L; ++l) {
            for (int m = 0; m <= l; ++m) {
                int idx = shIdx(l, m);
                double Plm = P[idx];
                double Nlm = normFactors[idx];
                double Ybase = Nlm * Plm;

                // cos branch (m ≥ 0): Y_l^m = N_l^m * P_l^m * cos(mφ)
                double Ycos = Ybase * cosMPhi[m];
                sumWU_cosM[idx] += wu * Ycos;
                sumW_cosM[idx] += w * Ycos;

                // sin branch (m > 0): Y_l^{-m} = N_l^m * P_l^m * sin(mφ)
                if (m > 0) {
                    double Ysin = Ybase * sinMPhi[m];
                    sumWU_sinM[idx] += wu * Ysin;
                    sumW_sinM[idx] += w * Ysin;
                }
            }
        }
        weightSum += w;
        sampleCount += 1;
    }

    void finalizeCoefficients(bool useSN, bool useCV)
    {
        if (sampleCount == 0)
            return;
        double denom = useSN ? weightSum : static_cast<double>(sampleCount);

        // DC coefficient with 4π scaling: a_0^0 = 4π · Σ w·u·y_0^0 / denom = ∫ u·y_0^0 dΩ
        ScalarType a00 = sumWU_cosM[shIdx(0, 0)] / denom * 4.0 * PI;
        aCoeffs[shIdx(0, 0)] = a00;
        bCoeffs[shIdx(0, 0)] = ScalarType(0.0);

        // Mean boundary value for CV: ū = a_0^0 · y_0^0 = a_0^0 / √(4π)
        ScalarType uMean = a00 * normFactors[shIdx(0, 0)];

        for (int l = 1; l <= L; ++l) {
            for (int m = 0; m <= l; ++m) {
                int idx = shIdx(l, m);
                ScalarType coeff = sumWU_cosM[idx] / denom * 4.0 * PI;
                if (useCV) {
                    double basisMean = sumW_cosM[idx] / denom * 4.0 * PI;
                    aCoeffs[idx] = coeff - uMean * basisMean;
                } else {
                    aCoeffs[idx] = coeff;
                }

                if (m > 0) {
                    ScalarType coeffSin = sumWU_sinM[idx] / denom * 4.0 * PI;
                    if (useCV) {
                        double basisMeanSin = sumW_sinM[idx] / denom * 4.0 * PI;
                        bCoeffs[idx] = coeffSin - uMean * basisMeanSin;
                    } else {
                        bCoeffs[idx] = coeffSin;
                    }
                }
            }
        }

        if (srcCount > 0)
            for (auto& s : srcAccum)
                s /= srcCount;
    }

    ScalarType evaluate(const Vector<3>& p, double kappa) const
    {
        Vector<3> offset = p - center;
        double d = offset.norm();
        if (d < EPSILON) {
            // Only l=0 contributes at center: u = a_0^0 · R_0(0) · y_0^0
            double R0 = radialBasis<3>(0.0, radius, 0, kappa);
            return aCoeffs[shIdx(0, 0)] * R0 * normFactors[shIdx(0, 0)];
        }

        double cosTheta = offset[2] / d;
        double phi = std::atan2(offset[1], offset[0]);

        // cos(mφ) and sin(mφ) via Chebyshev recurrence
        std::vector<double> cosMPhi(L + 1), sinMPhi(L + 1);
        cosMPhi[0] = 1.0;
        sinMPhi[0] = 0.0;
        if (L >= 1) {
            double cosPhi = std::cos(phi);
            double sinPhi = std::sin(phi);
            cosMPhi[1] = cosPhi;
            sinMPhi[1] = sinPhi;
            for (int m = 2; m <= L; ++m) {
                cosMPhi[m] = 2.0 * cosPhi * cosMPhi[m - 1] - cosMPhi[m - 2];
                sinMPhi[m] = 2.0 * cosPhi * sinMPhi[m - 1] - sinMPhi[m - 2];
            }
        }

        // Associated Legendre P_l^m(cosθ)
        std::vector<double> P(nCoeff, 0.0);
        double sx = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
        P[shIdx(0, 0)] = 1.0;

        for (int m = 1; m <= L; ++m) {
            P[shIdx(m, m)] = (2.0 * m - 1.0) * sx * P[shIdx(m - 1, m - 1)];
        }

        for (int m = 0; m <= L; ++m) {
            if (m + 1 <= L) {
                P[shIdx(m + 1, m)] = (2.0 * m + 1.0) * cosTheta * P[shIdx(m, m)];
            }
            for (int l = m + 2; l <= L; ++l) {
                double term1 = (2.0 * l - 1.0) * cosTheta * P[shIdx(l - 1, m)];
                double term2 = (l + m - 1.0) * P[shIdx(l - 2, m)];
                P[shIdx(l, m)] = (term1 - term2) / (l - m);
            }
        }

        ScalarType val(0.0);
        for (int l = 0; l <= L; ++l) {
            double Rl = radialBasis<3>(d, radius, l, kappa);
            ScalarType inner(0.0);
            for (int m = 0; m <= l; ++m) {
                int idx = shIdx(l, m);
                double Plm = P[idx];
                double Nlm = normFactors[idx];
                double Ybase = Nlm * Plm;

                // cos branch: a_l^m * Y_l^m
                inner += aCoeffs[idx] * (Ybase * cosMPhi[m]);

                // sin branch: b_l^m * Y_l^{-m}
                if (m > 0)
                    inner += bCoeffs[idx] * (Ybase * sinMPhi[m]);
            }
            val += Rl * inner;
        }
        if (val.isNaN()) {
            spdlog::error("Probe evaluation returned NaN at p={}, center={}, radius={}", p, center, radius);
            return ScalarType::NaN();
        }
        return val;
    }
};

WOS_NAMESPACE_CLOSE_SCOPE
