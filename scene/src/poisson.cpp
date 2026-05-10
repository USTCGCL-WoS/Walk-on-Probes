#include "scene/poisson.hpp"
#include "core/scalar.hpp"

WOS_NAMESPACE_OPEN_SCOPE

// 显式实例化
template class WOS_SCENE_API PoissonScene<Scalar<1>, 2>;
template class WOS_SCENE_API PoissonScene<Scalar<1>, 3>;
template class WOS_SCENE_API PoissonScene<Scalar<3>, 2>;
template class WOS_SCENE_API PoissonScene<Scalar<3>, 3>;
WOS_NAMESPACE_CLOSE_SCOPE
