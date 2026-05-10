#pragma once

#include "core/core_api.h"
#include <Eigen/Dense>
#include <spdlog/fmt/fmt.h>

WOS_NAMESPACE_OPEN_SCOPE

template<int DIM>
using Vector = Eigen::Matrix<double, DIM, 1>;

WOS_NAMESPACE_CLOSE_SCOPE

namespace fmt {

template<int DIM>
struct formatter<Eigen::Matrix<double, DIM, 1>>
{
    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    template<typename FormatContext>
    auto format(const Eigen::Matrix<double, DIM, 1>& v, FormatContext& ctx) const
    {
        auto out = ctx.out();
        *out++ = '[';
        for (int i = 0; i < DIM; ++i) {
            if (i > 0) {
                *out++ = ',';
                *out++ = ' ';
            }
            out = format_to(out, "{}", v[i]);
        }
        *out++ = ']';
        return out;
    }
};

} // namespace fmt