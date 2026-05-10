#pragma once
#include <numbers>

#include "core/core_api.h"

WOS_NAMESPACE_OPEN_SCOPE

inline constexpr double E = std::numbers::e;
inline constexpr double LOG2E = std::numbers::log2e;
inline constexpr double LOG10E = std::numbers::log10e;
inline constexpr double PI = std::numbers::pi;
inline constexpr double INV_PI = std::numbers::inv_pi;
inline constexpr double INV_SQRTPI = std::numbers::inv_sqrtpi;
inline constexpr double LN2 = std::numbers::ln2;
inline constexpr double LN10 = std::numbers::ln10;
inline constexpr double SQRT2 = std::numbers::sqrt2;
inline constexpr double SQRT3 = std::numbers::sqrt3;
inline constexpr double INV_SQRT3 = std::numbers::inv_sqrt3;
inline constexpr double EGAMMA = std::numbers::egamma;
inline constexpr double PHI = std::numbers::phi;

inline constexpr double EPSILON = 1e-4;
inline constexpr double EPSILON_SQ = EPSILON * EPSILON;

WOS_NAMESPACE_CLOSE_SCOPE