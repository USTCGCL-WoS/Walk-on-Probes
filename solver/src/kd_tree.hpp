#pragma once

#include "core/vector.hpp"
#include "solver/solver_api.h"

#include <nanoflann.hpp>

WOS_NAMESPACE_OPEN_SCOPE

template<int DIM>
struct PointCloudAdaptor
{
    const std::vector<Vector<DIM>>& points;

    explicit PointCloudAdaptor(const std::vector<Vector<DIM>>& pts)
      : points(pts)
    {
    }

    size_t kdtree_get_point_count() const { return points.size(); }

    double kdtree_get_pt(const size_t idx, const size_t dim) const { return points[idx][static_cast<int>(dim)]; }

    template<class BBOX>
    bool kdtree_get_bbox(BBOX&) const
    {
        return false;
    }
};

WOS_NAMESPACE_CLOSE_SCOPE