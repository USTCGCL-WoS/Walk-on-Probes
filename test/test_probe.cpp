#include "solver/probe.hpp"

#include "pcg_random.hpp"

#include <cmath>
#include <random>

#include <gtest/gtest.h>

using namespace WOS;

namespace {

// ---- typed test parameters ------------------------------------------------

template<int D>
struct Dim
{
    static constexpr int value = D;
};
using Dims = ::testing::Types<Dim<2>, Dim<3>>;

constexpr int kExpansionOrder2D = 10;
constexpr int kExpansionOrder3D = 4;
constexpr int kSampleCount = 20000;
constexpr double kR = 1.0;

double
tolerance()
{
    return 8.0 / std::sqrt(static_cast<double>(kSampleCount));
}

// =========================================================================
// Analytical test functions.
// Probe is at origin, radius kR.
// Test checks: probe.evaluate(p, κ) ≈ exactFunc(p)  for |p| < kR.
// =========================================================================

// ---- constant (κ = 0 only) -----------------------------------------------

template<int DIM>
double
constantExact(Vector<DIM> /*x*/)
{
    return 3.0;
}

// ---- linear harmonic (κ = 0 only) ----------------------------------------
// 2D: u = x + 2y     3D: u = x + 2y - 3z

template<int DIM>
double
linearExact(Vector<DIM> x)
{
    double val = x[0] + 2.0 * x[1];
    if constexpr (DIM == 3)
        val -= 3.0 * x[2];
    return val;
}

// ---- quadratic harmonic (κ = 0 only) -------------------------------------
// 2D: u = x² - y²     3D: u = 2z² - x² - y²

template<int DIM>
double
quadraticExact(Vector<DIM> x)
{
    if constexpr (DIM == 2)
        return x[0] * x[0] - x[1] * x[1];
    else
        return 2.0 * x[2] * x[2] - x[0] * x[0] - x[1] * x[1];
}

// ---- exponential screened (κ > 0) ----------------------------------------
// Δ(e^{κx}) - κ²(e^{κx}) = 0  in any dimension.

template<int DIM>
double
exponentialExact(Vector<DIM> x, double kappa)
{
    return std::exp(kappa * x[0]);
}

// =========================================================================
// Test harness
// =========================================================================

template<int DIM>
void
runTest(const std::function<double(Vector<DIM>)>& exactFunc,
        double kappa,
        const std::string& label)
{
    constexpr int L = (DIM == 2) ? kExpansionOrder2D : kExpansionOrder3D;

    Probe<Scalar<1>, DIM> probe;
    probe.center = Vector<DIM>::Zero();
    probe.radius = kR;
    probe.initExpansion(L);

    pcg64 rng(42);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    double pdf = (DIM == 2) ? 0.5 * INV_PI : 0.25 * INV_PI;

    for (int i = 0; i < kSampleCount; ++i) {
        Vector<DIM> dir;
        if constexpr (DIM == 2) {
            double theta = 2.0 * PI * u01(rng);
            dir << std::cos(theta), std::sin(theta);
        } else {
            double cosTheta = 1.0 - 2.0 * u01(rng);
            double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
            double phi = 2.0 * PI * u01(rng);
            dir << sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta;
        }
        Vector<DIM> samplePos = dir * kR;

        probe.addSampleUnlocked(samplePos, Scalar<1>(exactFunc(samplePos)), pdf);
    }

    probe.finalizeCoefficients(true, true);

    // Center
    {
        Vector<DIM> p = Vector<DIM>::Zero();
        Scalar<1> result = probe.evaluate(p, kappa);
        EXPECT_NEAR(result.value(), exactFunc(p), tolerance())
          << label << " center";
    }

    // Random interior points  (uniform in ball, capped at 0.9·R)
    pcg64 rng2(99);
    for (int i = 0; i < 8; ++i) {
        Vector<DIM> dir;
        double r;
        if constexpr (DIM == 2) {
            double theta = 2.0 * PI * u01(rng2);
            dir << std::cos(theta), std::sin(theta);
            r = std::sqrt(u01(rng2)) * kR * 0.9;
        } else {
            double cosTheta = 1.0 - 2.0 * u01(rng2);
            double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
            double phi = 2.0 * PI * u01(rng2);
            dir << sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta;
            r = std::cbrt(u01(rng2)) * kR * 0.9;
        }

        Vector<DIM> p = dir * r;
        Scalar<1> result = probe.evaluate(p, kappa);
        EXPECT_NEAR(result.value(), exactFunc(p), tolerance())
          << label << "  r=" << r;
    }
}

// Overload for exponential (needs kappa)
template<int DIM>
void
runScreenedTest(double kappa)
{
    auto f = [kappa](Vector<DIM> x) { return exponentialExact<DIM>(x, kappa); };
    runTest<DIM>(f, kappa, "Exponential κ=" + std::to_string(kappa));
}

} // namespace

// ---- typed tests ----------------------------------------------------------

template<typename T>
class ProbeConstantTest : public ::testing::Test
{};
TYPED_TEST_SUITE(ProbeConstantTest, Dims);
TYPED_TEST(ProbeConstantTest, Laplace)
{
    runTest<TypeParam::value>(constantExact<TypeParam::value>, 0.0, "Constant");
}

template<typename T>
class ProbeLinearTest : public ::testing::Test
{};
TYPED_TEST_SUITE(ProbeLinearTest, Dims);
TYPED_TEST(ProbeLinearTest, Laplace)
{
    runTest<TypeParam::value>(linearExact<TypeParam::value>, 0.0, "Linear");
}

template<typename T>
class ProbeQuadraticTest : public ::testing::Test
{};
TYPED_TEST_SUITE(ProbeQuadraticTest, Dims);
TYPED_TEST(ProbeQuadraticTest, Laplace)
{
    runTest<TypeParam::value>(quadraticExact<TypeParam::value>, 0.0, "Quadratic");
}

template<typename T>
class ProbeScreenedTest : public ::testing::Test
{};
TYPED_TEST_SUITE(ProbeScreenedTest, Dims);
TYPED_TEST(ProbeScreenedTest, Exponential)
{
    runScreenedTest<TypeParam::value>(0.5);
}
