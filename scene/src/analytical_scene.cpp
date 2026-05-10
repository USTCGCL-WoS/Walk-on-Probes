#include "scene/analytical_scene.hpp"

#include "core/math_defs.hpp"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
AnalyticalScene<ScalarType, DIM>::AnalyticalScene(const std::string& sceneDir)
  : ObjScene<ScalarType, DIM>(sceneDir)
{
}

template<typename ScalarType, int DIM>
ScalarType
AnalyticalScene<ScalarType, DIM>::exactSolution(const VectorType& p) const
{
    ScalarType v(1.0);
    for (int d = 0; d < DIM; ++d) {
        v *= std::sin(PI * p[d]);
    }
    return (v + 1.0) * 0.5;
}

template<typename ScalarType, int DIM>
ScalarType
AnalyticalScene<ScalarType, DIM>::exactGradientDotNormal(const VectorType& p, const VectorType& n) const
{
    ScalarType result(0.0);
    for (int d = 0; d < DIM; ++d) {
        ScalarType grad_d(1.0);
        for (int k = 0; k < DIM; ++k) {
            if (k == d)
                grad_d *= PI * std::cos(PI * p[k]);
            else
                grad_d *= std::sin(PI * p[k]);
        }
        result += grad_d * n[d];
    }
    return result * 0.5;
}

template<typename ScalarType, int DIM>
ScalarType
AnalyticalScene<ScalarType, DIM>::source(const VectorType& p) const
{
    double kappa = this->isScreened() ? this->getScreenedKappa() : 0.0;
    double Dpi2 = DIM * PI * PI;
    return exactSolution(p) * (Dpi2 + kappa) - ScalarType(Dpi2 * 0.5);
}

template<typename ScalarType, int DIM>
BoundaryCondition<ScalarType>
AnalyticalScene<ScalarType, DIM>::boundaryCondition(const fcpw::Interaction<DIM>& interaction,
                                                    const VectorType& queryPoint,
                                                    BoundaryType type) const
{
    VectorType p = interaction.p.template cast<double>();

    switch (type) {
        case BoundaryType::Dirichlet:
            return BoundaryCondition<ScalarType>::Dirichlet(exactSolution(p));
        case BoundaryType::Neumann: {
            VectorType n = interaction.n.template cast<double>();
            ScalarType val = exactGradientDotNormal(p, n);
            if ((queryPoint.template cast<float>() - interaction.p).dot(interaction.n) > 0)
                val = -val;
            return BoundaryCondition<ScalarType>::Neumann(val);
        }
        default:
            return BoundaryCondition<ScalarType>();
    }
}

template class WOS_SCENE_API AnalyticalScene<Scalar<1>, 2>;
template class WOS_SCENE_API AnalyticalScene<Scalar<1>, 3>;
template class WOS_SCENE_API AnalyticalScene<Scalar<3>, 2>;
template class WOS_SCENE_API AnalyticalScene<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
