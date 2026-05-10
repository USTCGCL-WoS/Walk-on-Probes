#pragma once

#include <Eigen/Core>
#include <memory>

#include "scene/scene_api.h"

#include <utility>
#include <vector>

#include "core/math_defs.hpp"
#include "core/scalar.hpp"
#include "core/vector.hpp"
#include "fcpw/core/interaction.h"
#include "fcpw/fcpw.h"
#include "pcg_random.hpp"
#include <random>

WOS_NAMESPACE_OPEN_SCOPE

enum class BoundaryType
{
    Dirichlet,
    Neumann,
    Robin
};

template<typename ScalarType>
struct BoundaryCondition
{
    BoundaryType type;
    ScalarType a; // Dirichlet: boundary value; Neumann: unused
    ScalarType b; // Neumann: normal derivative; Dirichlet: unused

    BoundaryCondition() noexcept
      : type(BoundaryType::Dirichlet)
      , a(ScalarType(0.0))
      , b(ScalarType(0.0))
    {
    }

    BoundaryCondition(BoundaryType t, const ScalarType& a_val, const ScalarType& b_val) noexcept
      : type(t)
      , a(a_val)
      , b(b_val)
    {
    }

    [[nodiscard]] static inline BoundaryCondition Dirichlet(const ScalarType& value) noexcept
    {
        return BoundaryCondition(BoundaryType::Dirichlet, value, ScalarType(0.0));
    }

    [[nodiscard]] static inline BoundaryCondition Neumann(const ScalarType& value) noexcept
    {
        return BoundaryCondition(BoundaryType::Neumann, ScalarType(0.0), value);
    }
};

/// Screened Poisson equation: Δu - κu = f in Ω with boundary conditions.
template<typename ScalarType, int DIM>
class PoissonScene
{
  public:
    PoissonScene() {}
    virtual ~PoissonScene() = default;

    virtual ScalarType source(const Vector<DIM>& p) const = 0;

    virtual BoundaryCondition<ScalarType> boundaryCondition(const fcpw::Interaction<DIM>& interaction,
                                                            const Vector<DIM>& queryPoint,
                                                            BoundaryType type) const = 0;

    bool isScreened() const { return is_screened_; }

    double getScreenedKappa() const { return screened_kappa_; }

    bool findClosestBoundary(const Vector<DIM>& p, fcpw::Interaction<DIM>& interaction) const
    {
        if (!full_scene_)
            return false;
        return full_scene_->findClosestPoint(p.template cast<float>(), interaction, fcpw::maxFloat, true);
    }

    bool findClosestDirichletBoundary(const Vector<DIM>& p, fcpw::Interaction<DIM>& interaction) const
    {
        if (!dirichlet_scene_)
            return false;
        return dirichlet_scene_->findClosestPoint(p.template cast<float>(), interaction, fcpw::maxFloat, true);
    }

    bool findClosestNeumannBoundary(const Vector<DIM>& p, fcpw::Interaction<DIM>& interaction) const
    {
        if (!neumann_scene_)
            return false;
        return neumann_scene_->findClosestPoint(p.template cast<float>(), interaction, fcpw::maxFloat, true);
    }

    bool findClosestNeumannSilhouette(const Vector<DIM>& p,
                                      fcpw::Interaction<DIM>& interaction,
                                      bool flipNormalOrientation,
                                      float rMinSquared,
                                      float rMaxSquared,
                                      float precision) const
    {
        if (!neumann_scene_)
            return false;
        return neumann_scene_->findClosestSilhouettePoint(
          p.template cast<float>(), interaction, flipNormalOrientation, rMinSquared, rMaxSquared, precision, true);
    }

    bool rayIntersectNeumann(const Vector<DIM>& origin,
                             const Vector<DIM>& dir,
                             float max_dist,
                             fcpw::Interaction<DIM>& interaction) const
    {
        if (!neumann_scene_)
            return false;
        fcpw::Ray<DIM> ray(origin.template cast<float>(), dir.template cast<float>(), max_dist);
        return neumann_scene_->intersect(ray, interaction, false);
    }

    bool sampleNeumannBoundary(const Vector<DIM>& center,
                               float radius,
                               pcg64& rng,
                               std::uniform_real_distribution<double>& dist01,
                               fcpw::Interaction<DIM>& interaction) const
    {
        if (!neumann_scene_)
            return false;

        auto weight_fn = [](float r2) -> float {
            float r = std::sqrt(r2);
            r = std::max(r, static_cast<float>(EPSILON));
            return 0.25f / (PI * r);
        };

        fcpw::BoundingSphere<DIM> sphere(center.template cast<float>(), radius * radius);
        fcpw::Vector<DIM> randNums;
        for (int d = 0; d < DIM; ++d) {
            randNums[d] = dist01(rng);
        }

        int hits = this->neumann_scene_->intersect(sphere, interaction, randNums, weight_fn);
        return hits > 0;
    }

    double distanceToInterestBoundary(const Vector<DIM>& p, fcpw::Interaction<DIM>& interaction) const
    {
        if (!isInsideROI(p))
            return -1.0;

        if (has_interest_) {
            interest_scene_->findClosestPoint(p.template cast<float>(), interaction, fcpw::maxFloat, false);
            return interaction.d;
        }

        double min_dist = std::numeric_limits<double>::max();
        for (int d = 0; d < DIM; ++d) {
            double dist_min = p[d] - renderMin_[d];
            double dist_max = renderMax_[d] - p[d];
            min_dist = std::min({ min_dist, dist_min, dist_max });
        }
        return min_dist;
    }

    bool checkLineOfSight(const Vector<DIM>& p1, const Vector<DIM>& p2) const
    {
        return this->neumann_scene_->hasLineOfSight(p1.template cast<float>(), p2.template cast<float>());
    }

    virtual bool isInsideROI(const Vector<DIM>& p) const
    {
        bool inside = true;
        if (interior_only_ && full_scene_) {
            inside = full_scene_->contains(p.template cast<float>());
        }
        if (has_interest_) {
            inside = inside && interest_scene_->contains(p.template cast<float>());
        }
        // for (int d = 0; d < DIM; ++d) {
        //     if (p[d] < renderMin_[d] || p[d] > renderMax_[d])
        //         return false;
        // }
        return inside;
    }

    bool isInteriorOnly() const { return interior_only_; }

    virtual bool isInsideObject(const Vector<DIM>& p) const
    {
        if (full_scene_) {
            return full_scene_->contains(p.template cast<float>());
        }
        return true;
    }

    Vector<DIM> getRenderMin() const { return renderMin_; }
    Vector<DIM> getRenderMax() const { return renderMax_; }

    virtual std::pair<std::vector<Vector<DIM>>, std::vector<std::size_t>> generateTargetPoints(int width,
                                                                                               int height) const
    {
        return {};
    }

  protected:
    std::unique_ptr<fcpw::Scene<DIM>> dirichlet_scene_{ nullptr };
    std::unique_ptr<fcpw::Scene<DIM>> neumann_scene_{ nullptr };
    std::unique_ptr<fcpw::Scene<DIM>> interest_scene_{ nullptr };
    std::unique_ptr<fcpw::Scene<DIM>> full_scene_{ nullptr };

    bool interior_only_ = false;
    bool has_interest_ = false;
    bool is_screened_ = false;

    double screened_kappa_ = 0.0;

    Vector<DIM> renderMin_{ Vector<DIM>() };
    Vector<DIM> renderMax_{ Vector<DIM>() };
};

extern template class PoissonScene<Scalar<1>, 2>;
extern template class PoissonScene<Scalar<1>, 3>;
extern template class PoissonScene<Scalar<3>, 2>;
extern template class PoissonScene<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
