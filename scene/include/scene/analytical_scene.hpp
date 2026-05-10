#pragma once

#include "scene/obj_scene.hpp"
#include "scene/scene_api.h"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
class AnalyticalScene : public ObjScene<ScalarType, DIM>
{
  public:
    using VectorType = Vector<DIM>;

    explicit AnalyticalScene(const std::string& sceneDir);

    ScalarType source(const VectorType& p) const override;

    BoundaryCondition<ScalarType> boundaryCondition(const fcpw::Interaction<DIM>& interaction,
                                                    const VectorType& queryPoint,
                                                    BoundaryType type) const override;

    ScalarType exactSolution(const VectorType& p) const;
    ScalarType exactGradientDotNormal(const VectorType& p, const VectorType& n) const;
};

extern template class AnalyticalScene<Scalar<1>, 2>;
extern template class AnalyticalScene<Scalar<1>, 3>;
extern template class AnalyticalScene<Scalar<3>, 2>;
extern template class AnalyticalScene<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
