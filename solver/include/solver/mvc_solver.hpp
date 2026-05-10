#pragma once

#include <memory>
#include <vector>

#include "core/vector.hpp"
#include "solver/caching_solver.hpp"
#include "solver/solver_api.h"
#include "solver/wost_solver.hpp"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
class MVCSolver : public CachingSolver<ScalarType, DIM>
{
  public:
    using CachingSolver<ScalarType, DIM>::CachingSolver;

    void buildCache() override;

    void configure(const nlohmann::json& j) override;

    ScalarType solve(const Vector<DIM>& targetPoint) override;

    void solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results) override;

    // Configuration setters
    void setNumIterations(int n) { numIterations_ = n; }
    void setNumCachePoints(int n) { numCachePoints_ = n; }
    void setWoStWalksPerPixel(int wpp) { wostWpp_ = wpp; }
    void setMaxSteps(int n) { maxSteps_ = n; }
    void setEpsilon(double eps) { epsilon_ = eps; }
    void setRMin(double r) { rMin_ = r; }
    void setSourceRatio(double r) { sourceRatio_ = r; }
    void setMinSourceSamples(int n) { minSourceSamples_ = n; }
    void setMu(double m) { mu_ = m; }
    void setLambda(double l) { lambda_ = l; }
    void setAdaptQ(double q) { adaptQ_ = q; }

  private:
    struct CachePoint
    {
        Vector<DIM> position;
        ScalarType solution;
        double importanceWeight = 1.0;
    };

    // Per-round helpers
    void generateCachePoints(std::vector<CachePoint>& out);

    // Query point state
    std::vector<Vector<DIM>> targetPoints_;
    std::vector<ScalarType> accumulatedSolutions_;
    std::vector<int> accumulatedCounts_;
    bool cacheBuilt_ = false;

    // Delegate solver
    std::unique_ptr<WoStSolver<ScalarType, DIM>> woStSolver_;

    // Configuration
    int numIterations_ = 4;
    int numCachePoints_ = 256;
    int wostWpp_ = 1;
    int maxSteps_ = 1024;
    double epsilon_ = EPSILON;
    double rMin_ = EPSILON;
    double sourceRatio_ = 1.0;
    int minSourceSamples_ = 1;
    double mu_ = 0.02;
    double lambda_ = 0.02;
    double adaptQ_ = 1.0;
};

extern template class MVCSolver<Scalar<1>, 2>;
extern template class MVCSolver<Scalar<1>, 3>;
extern template class MVCSolver<Scalar<3>, 2>;
extern template class MVCSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
