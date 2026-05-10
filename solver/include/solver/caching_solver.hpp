#pragma once

#include "solver/solver.hpp"
#include "solver/solver_api.h"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
class CachingSolver : public Solver<ScalarType, DIM>
{
  public:
    using Solver<ScalarType, DIM>::Solver;

    virtual void buildCache() = 0;
};

extern template class CachingSolver<Scalar<1>, 2>;
extern template class CachingSolver<Scalar<1>, 3>;
extern template class CachingSolver<Scalar<3>, 2>;
extern template class CachingSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
