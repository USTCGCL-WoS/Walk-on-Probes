#pragma once

#include "macros.h"
#include <Eigen/Dense>
#include <cmath>
#include <format>
#include <spdlog/fmt/fmt.h>

WOS_NAMESPACE_OPEN_SCOPE

template<int N>
class Scalar
{
    static_assert(N > 0, "Scalar<N> requires N > 0");

  public:
    static constexpr int channels = N;

    Scalar() = default;

    // Broadcast fill: Scalar<3>(1.0) -> [1, 1, 1]
    /* implicit */ Scalar(double v)
      : data(v)
    {
    }

    // 3-channel convenience constructor
    Scalar(double x, double y, double z)
        requires(N == 3)
      : data(x, y, z)
    {
    }

    // Element access
    double operator[](int i) const { return data[i]; }
    double& operator[](int i) { return data[i]; }

    // Raw Eigen access
    const auto& value() const { return data; }
    auto& value() { return data; }

    // ---- element-wise arithmetic ----
    Scalar operator+(const Scalar& o) const
    {
        Scalar r;
        r.data = data + o.data;
        return r;
    }
    Scalar operator-(const Scalar& o) const
    {
        Scalar r;
        r.data = data - o.data;
        return r;
    }
    Scalar operator*(const Scalar& o) const
    {
        Scalar r;
        r.data = data * o.data;
        return r;
    }
    Scalar operator/(const Scalar& o) const
    {
        Scalar r;
        r.data = data / o.data;
        return r;
    }

    Scalar& operator+=(const Scalar& o)
    {
        data += o.data;
        return *this;
    }
    Scalar& operator-=(const Scalar& o)
    {
        data -= o.data;
        return *this;
    }
    Scalar& operator*=(const Scalar& o)
    {
        data *= o.data;
        return *this;
    }
    Scalar& operator/=(const Scalar& o)
    {
        data /= o.data;
        return *this;
    }

    // ---- scalar (double) arithmetic ----
    Scalar operator+(double s) const
    {
        Scalar r;
        r.data = data + s;
        return r;
    }
    Scalar operator-(double s) const
    {
        Scalar r;
        r.data = data - s;
        return r;
    }
    Scalar operator*(double s) const
    {
        Scalar r;
        r.data = data * s;
        return r;
    }
    Scalar operator/(double s) const
    {
        Scalar r;
        r.data = data / s;
        return r;
    }

    Scalar& operator+=(double s)
    {
        data += s;
        return *this;
    }
    Scalar& operator-=(double s)
    {
        data -= s;
        return *this;
    }
    Scalar& operator*=(double s)
    {
        data *= s;
        return *this;
    }
    Scalar& operator/=(double s)
    {
        data /= s;
        return *this;
    }

    friend Scalar operator+(double s, const Scalar& v)
    {
        Scalar r;
        r.data = s + v.data;
        return r;
    }
    friend Scalar operator-(double s, const Scalar& v)
    {
        Scalar r;
        r.data = s - v.data;
        return r;
    }
    friend Scalar operator*(double s, const Scalar& v)
    {
        Scalar r;
        r.data = s * v.data;
        return r;
    }
    friend Scalar operator/(double s, const Scalar& v)
    {
        Scalar r;
        r.data = s / v.data;
        return r;
    }

    // Unary
    Scalar operator-() const
    {
        Scalar r;
        r.data = -data;
        return r;
    }
    Scalar operator+() const { return *this; }

    // Comparisons
    bool operator==(const Scalar& o) const { return (data == o.data).all(); }
    bool operator!=(const Scalar& o) const { return !(*this == o); }

    // ---- component-wise math ----
    Scalar cwiseMin(const Scalar& o) const
    {
        Scalar r;
        r.data = data.min(o.data);
        return r;
    }
    Scalar cwiseMax(const Scalar& o) const
    {
        Scalar r;
        r.data = data.max(o.data);
        return r;
    }
    Scalar cwiseAbs() const
    {
        Scalar r;
        r.data = data.abs();
        return r;
    }
    Scalar cwiseSqrt() const
    {
        Scalar r;
        r.data = data.sqrt();
        return r;
    }

    // ---- reductions ----
    double sum() const { return data.sum(); }
    double prod() const { return data.prod(); }
    double mean() const { return data.mean(); }
    double minCoeff() const { return data.minCoeff(); }
    double maxCoeff() const { return data.maxCoeff(); }

    // ---- NaN Utility ----
    static Scalar NaN()
    {
        Scalar r;
        r.data.setConstant(std::numeric_limits<double>::quiet_NaN());
        return r;
    }

    bool isNaN() const { return data.hasNaN(); }

  private:
    Eigen::Array<double, N, 1> data = Eigen::Array<double, N, 1>::Zero();
};

// ============================================================
// N=1 specialization — backed by plain double
// ============================================================
template<>
class Scalar<1>
{
  public:
    static constexpr int channels = 1;

    Scalar() = default;
    /* implicit */ Scalar(double v)
      : data(v)
    {
    }
    /* implicit */ operator double() const { return data; }

    double operator[](int) const { return data; }
    double& operator[](int) { return data; }

    double value() const { return data; }
    double& value() { return data; }

    Scalar operator+(Scalar o) const { return Scalar(data + o.data); }
    Scalar operator-(Scalar o) const { return Scalar(data - o.data); }
    Scalar operator*(Scalar o) const { return Scalar(data * o.data); }
    Scalar operator/(Scalar o) const { return Scalar(data / o.data); }

    Scalar& operator+=(Scalar o)
    {
        data += o.data;
        return *this;
    }
    Scalar& operator-=(Scalar o)
    {
        data -= o.data;
        return *this;
    }
    Scalar& operator*=(Scalar o)
    {
        data *= o.data;
        return *this;
    }
    Scalar& operator/=(Scalar o)
    {
        data /= o.data;
        return *this;
    }

    Scalar operator+(double s) const { return Scalar(data + s); }
    Scalar operator-(double s) const { return Scalar(data - s); }
    Scalar operator*(double s) const { return Scalar(data * s); }
    Scalar operator/(double s) const { return Scalar(data / s); }

    Scalar& operator+=(double s)
    {
        data += s;
        return *this;
    }
    Scalar& operator-=(double s)
    {
        data -= s;
        return *this;
    }
    Scalar& operator*=(double s)
    {
        data *= s;
        return *this;
    }
    Scalar& operator/=(double s)
    {
        data /= s;
        return *this;
    }

    friend Scalar operator+(double s, Scalar v) { return Scalar(s + v.data); }
    friend Scalar operator-(double s, Scalar v) { return Scalar(s - v.data); }
    friend Scalar operator*(double s, Scalar v) { return Scalar(s * v.data); }
    friend Scalar operator/(double s, Scalar v) { return Scalar(s / v.data); }

    Scalar operator-() const { return Scalar(-data); }
    Scalar operator+() const { return *this; }

    bool operator==(Scalar o) const { return data == o.data; }
    bool operator!=(Scalar o) const { return data != o.data; }

    Scalar cwiseMin(Scalar o) const { return Scalar(std::fmin(data, o.data)); }
    Scalar cwiseMax(Scalar o) const { return Scalar(std::fmax(data, o.data)); }
    Scalar cwiseAbs() const { return Scalar(std::abs(data)); }
    Scalar cwiseSqrt() const { return Scalar(std::sqrt(data)); }

    double sum() const { return data; }
    double prod() const { return data; }
    double mean() const { return data; }
    double minCoeff() const { return data; }
    double maxCoeff() const { return data; }

    static Scalar NaN() { return Scalar(std::numeric_limits<double>::quiet_NaN()); }

    bool isNaN() const { return std::isnan(data); }

  private:
    double data = 0.0;
};

// Type aliases
using Scalar1d = Scalar<1>;
using Scalar3d = Scalar<3>;

// ---- free functions ----
template<int N>
Scalar<N>
min(const Scalar<N>& a, const Scalar<N>& b)
{
    return a.cwiseMin(b);
}

template<int N>
Scalar<N>
max(const Scalar<N>& a, const Scalar<N>& b)
{
    return a.cwiseMax(b);
}

template<int N>
Scalar<N>
abs(const Scalar<N>& a)
{
    return a.cwiseAbs();
}

template<int N>
Scalar<N>
sqrt(const Scalar<N>& a)
{
    return a.cwiseSqrt();
}

template<int N>
double
dot(const Scalar<N>& a, const Scalar<N>& b)
{
    double sum = 0;
    for (int i = 0; i < N; ++i)
        sum += a[i] * b[i];
    return sum;
}

WOS_NAMESPACE_CLOSE_SCOPE

namespace fmt {

template<int N>
struct formatter<WOS::Scalar<N>>
{
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    template<typename FormatContext>
    auto format(const WOS::Scalar<N>& s, FormatContext& ctx) const
    {
        auto out = ctx.out();
        *out++ = '[';
        for (int i = 0; i < N; ++i) {
            if (i > 0) {
                *out++ = ',';
                *out++ = ' ';
            }
            out = format_to(out, "{}", s[i]);
        }
        *out++ = ']';
        return out;
    }
};

} // namespace fmt