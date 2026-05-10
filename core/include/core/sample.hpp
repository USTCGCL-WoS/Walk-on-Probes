#pragma once

#include "core/vector.hpp"

WOS_NAMESPACE_OPEN_SCOPE

// Result of a sampling operation: the sampled point and its probability density
template<int DIM>
struct Sample
{
    Vector<DIM> point;
    double pdf;
};

template<int DIM>
struct SampleLength
{
    double length;
    double pdf;
};

WOS_NAMESPACE_CLOSE_SCOPE
