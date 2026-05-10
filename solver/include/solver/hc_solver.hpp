#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include "core/vector.hpp"
#include "solver/caching_solver.hpp"
#include "solver/probe.hpp"
#include "solver/probe_set.hpp"
#include "solver/solver_api.h"
#include "solver/wost_solver.hpp"

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
class HCSolver : public CachingSolver<ScalarType, DIM>
{
  public:
    using CachingSolver<ScalarType, DIM>::CachingSolver;

    void buildCache() override;

    void configure(const nlohmann::json& j) override;

    ScalarType solve(const Vector<DIM>& targetPoint) override;

    void solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results) override;

    // Configuration setters
    void setTargetPoints(const std::vector<Vector<DIM>>& pts) { targetPoints_ = pts; }
    void setNumIterations(int n) { numIterations_ = n; }
    void setFourierModes(int L) { fourierModes_ = L; }
    void setSHModes(int L) { shModes_ = L; }
    void setBoundarySampleFactor(double lambda) { lambda_ = lambda; }
    void setSourceSampleFactor(double mu) { mu_ = mu; }
    void setMinBoundarySamples(int n) { minBoundarySamples_ = n; }
    void setMinSourceSamples(int n) { minSourceSamples_ = n; }
    void setAlphaRec(double a) { alphaRec_ = a; }
    void setWMin(double w) { wMin_ = w; }
    void setWoStWalksPerPixel(int wpp) { wostWpp_ = wpp; }
    void setMaxSteps(int n) { maxSteps_ = n; }
    void setEpsilon(double eps) { epsilon_ = eps; }
    void setRMin(double r) { rMin_ = r; }
    void setUseSN(bool v) { useSN_ = v; }
    void setUseCV(bool v) { useCV_ = v; }

  private:
    // Internal helpers
    void runIteration(int iter);
    ScalarType solvePoint(const Vector<DIM>& p,
                          int pointIdx,
                          pcg64& rng,
                          std::uniform_real_distribution<double>& dist,
                          fcpw::Interaction<static_cast<size_t>(DIM)>& interaction);
    void sampleBoundarySamples(const Probe<ScalarType, DIM>& probe, int N, std::vector<Vector<DIM>>& boundarySamples);

    // Cached data
    ProbeSet<ScalarType, DIM> probes_;
    std::vector<Vector<DIM>> targetPoints_;
    bool cacheBuilt_ = false;
    std::unique_ptr<WoStSolver<ScalarType, DIM>> woStSolver_;
    std::vector<Vector<DIM>> allBoundaryPoints_;
    std::vector<ScalarType> allBoundaryResults_;
    std::vector<int> probeBoundaryCounts_;     // number of boundary samples per probe
    std::vector<size_t> probeBoundaryOffsets_; // prefix sum: probeBoundaryOffsets_[p] = start index in flat arrays

    // Configuration
    int numIterations_ = 1;
    int fourierModes_ = 10;
    int shModes_ = 4;
    double lambda_ = 100.0;
    double mu_ = 1000.0;
    int minBoundarySamples_ = 32;
    int minSourceSamples_ = 1;
    double alphaRec_ = 0.9;
    double wMin_ = 1.0;
    int wostWpp_ = 1;
    int maxSteps_ = 1024;
    double epsilon_ = EPSILON;
    double rMin_ = EPSILON;
    bool useSN_ = true;
    bool useCV_ = true;
};

extern template class HCSolver<Scalar<1>, 2>;
extern template class HCSolver<Scalar<1>, 3>;
extern template class HCSolver<Scalar<3>, 2>;
extern template class HCSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
