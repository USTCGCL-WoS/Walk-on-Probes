#include "core/math_defs.hpp"
#include "core/sampling_green.hpp"
#include "core/sampling_poisson_kernel.hpp"

#include "pcg_random.hpp"

#include <cmath>
#include <random>

#include <gtest/gtest.h>

using namespace WOS;

// ---------- typed test infrastructure ----------
template<int N>
struct Dim
{
    static constexpr int value = N;
};
using Dims = ::testing::Types<Dim<2>, Dim<3>>;

// ---------- helpers ----------
namespace {

constexpr double kR = 1.0;
constexpr int kNumSamples = 30000;

double
tolerance(double scale)
{
    return 4.0 * scale / std::sqrt(static_cast<double>(kNumSamples));
}

template<int DIM>
double
surfaceArea(double R)
{
    if constexpr (DIM == 2)
        return 2.0 * PI * R;
    return 4.0 * PI * R * R;
}

template<int DIM>
double
ballVolume(double R)
{
    if constexpr (DIM == 2)
        return PI * R * R;
    return (4.0 / 3.0) * PI * R * R * R;
}

} // namespace

// ============================================================
// Poisson kernel — at center: uniform on sphere → E[z_x]=0
// ============================================================
template<typename T>
class PoissonCenterTest : public ::testing::Test
{};

TYPED_TEST_SUITE(PoissonCenterTest, Dims);

TYPED_TEST(PoissonCenterTest, Moments)
{
    constexpr int D = TypeParam::value;
    pcg64 rng(42);
    std::uniform_real_distribution<double> u(0, 1);

    double sum_x = 0.0, sum_x2 = 0.0;
    for (int i = 0; i < kNumSamples; ++i) {
        auto [z, pdf] = samplePoissonKernelSphereAtCenter<D>(rng, u, kR);
        sum_x += z.x();
        sum_x2 += z.x() * z.x();
    }
    EXPECT_NEAR(sum_x / kNumSamples, 0.0, tolerance(kR));
    EXPECT_NEAR(sum_x2 / kNumSamples, kR * kR / D, tolerance(kR * kR));
}

// ============================================================
// Poisson kernel — off center: Poisson integral E[h(z)] = h(x)
// ============================================================
template<typename T>
class PoissonOffcenterTest : public ::testing::Test
{};

TYPED_TEST_SUITE(PoissonOffcenterTest, Dims);

TYPED_TEST(PoissonOffcenterTest, PoissonIntegral)
{
    constexpr int D = TypeParam::value;
    pcg64 rng(99);
    std::uniform_real_distribution<double> u(0, 1);

    Vector<D> x = Vector<D>::Zero();
    x.x() = 0.5 * kR;

    double sum_h1 = 0.0, sum_h2 = 0.0, sum_inv_pdf = 0.0;
    for (int i = 0; i < kNumSamples; ++i) {
        auto [z, pdf] = samplePoissonKernelSphereOffCenter<D>(rng, u, x, kR);
        sum_h1 += z.x();
        if constexpr (D == 2) {
            sum_h2 += z.x() * z.x() - z.y() * z.y();
        } else {
            sum_h2 += 2.0 * z.z() * z.z() - z.x() * z.x() - z.y() * z.y();
        }
        sum_inv_pdf += 1.0 / pdf;
    }

    // h₂(x) = x² in 2D, -x² in 3D (for the chosen harmonic polynomial)
    const double h2_x = (D == 2) ? (x.x() * x.x()) : (-x.x() * x.x());

    EXPECT_NEAR(sum_h1 / kNumSamples, x.x(), tolerance(kR));
    EXPECT_NEAR(sum_h2 / kNumSamples, h2_x, tolerance(kR * kR));
    EXPECT_NEAR(sum_inv_pdf / kNumSamples, surfaceArea<D>(kR), surfaceArea<D>(kR) * 0.05);
}

// ============================================================
// Green's function — at center: radial moments
// ============================================================
template<typename T>
class GreensCenterTest : public ::testing::Test
{};

TYPED_TEST_SUITE(GreensCenterTest, Dims);

TYPED_TEST(GreensCenterTest, RadialMoment)
{
    constexpr int D = TypeParam::value;
    pcg64 rng(123);
    std::uniform_real_distribution<double> u(0, 1);

    double sum_r = 0.0, sum_inv_pdf = 0.0;
    for (int i = 0; i < kNumSamples; ++i) {
        auto [y, pdf] = sampleGreensBallAtCenter<D>(rng, u, kR);
        sum_r += y.norm();
        sum_inv_pdf += 1.0 / pdf;
    }

    if constexpr (D == 2) {
        EXPECT_NEAR(sum_r / kNumSamples, 4.0 * kR / 9.0, tolerance(kR));
    } else {
        EXPECT_NEAR(sum_r / kNumSamples, 0.5 * kR, tolerance(kR));
    }
    EXPECT_NEAR(sum_inv_pdf / kNumSamples, ballVolume<D>(kR), ballVolume<D>(kR) * 0.05);
}

// ============================================================
// Green's function — off center: all points inside ball
// ============================================================
template<typename T>
class GreensOffcenterTest : public ::testing::Test
{};

TYPED_TEST_SUITE(GreensOffcenterTest, Dims);

TYPED_TEST(GreensOffcenterTest, InsideBall)
{
    constexpr int D = TypeParam::value;
    pcg64 rng(77);
    std::uniform_real_distribution<double> u(0, 1);

    Vector<D> x = Vector<D>::Zero();
    x.x() = 0.5 * kR;

    int inside_count = 0;
    double sum_inv_pdf = 0.0;
    for (int i = 0; i < kNumSamples; ++i) {
        auto [y, pdf] = sampleGreensBallOffCenter<D>(rng, u, x, kR);
        if (y.norm() <= kR)
            ++inside_count;
        sum_inv_pdf += 1.0 / pdf;
    }

    EXPECT_EQ(inside_count, kNumSamples);
    EXPECT_NEAR(sum_inv_pdf / kNumSamples, ballVolume<D>(kR), ballVolume<D>(kR) * 0.05);
}