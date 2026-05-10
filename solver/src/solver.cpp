#include "solver/solver.hpp"
#include "solver/caching_solver.hpp"

WOS_NAMESPACE_OPEN_SCOPE

// Solver explicit instantiations
template class WOS_SOLVER_API Solver<Scalar<1>, 2>;
template class WOS_SOLVER_API Solver<Scalar<1>, 3>;
template class WOS_SOLVER_API Solver<Scalar<3>, 2>;
template class WOS_SOLVER_API Solver<Scalar<3>, 3>;

// CachingSolver explicit instantiations
template class WOS_SOLVER_API CachingSolver<Scalar<1>, 2>;
template class WOS_SOLVER_API CachingSolver<Scalar<1>, 3>;
template class WOS_SOLVER_API CachingSolver<Scalar<3>, 2>;
template class WOS_SOLVER_API CachingSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
