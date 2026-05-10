#include "solver/wost_solver.hpp"
#include "walk.hpp"

#include <nlohmann/json.hpp>

WOS_NAMESPACE_OPEN_SCOPE

template<typename ScalarType, int DIM>
void
WoStSolver<ScalarType, DIM>::configure(const nlohmann::json& j)
{
    if (j.contains("wpp"))
        setWalksPerPixel(j["wpp"].get<int>());
    if (j.contains("max_steps"))
        setMaxSteps(j["max_steps"].get<int>());
    if (j.contains("epsilon"))
        setEpsilon(j["epsilon"].get<double>());
    if (j.contains("r_min"))
        setRMin(j["r_min"].get<double>());
}

template<typename ScalarType, int DIM>
ScalarType
WoStSolver<ScalarType, DIM>::solve(const Vector<DIM>& targetPoint)
{
    ScalarType result(0.0);
    fcpw::Interaction<DIM> interaction;
    for (int w = 0; w < wpp_; ++w) {
        auto walkResult =
          walkWoSt(*this->scene_, targetPoint, maxSteps_, epsilon_, rMin_, this->rng_, this->dist01_, interaction);
        if (walkResult.isNaN()) {
            --w; // retry this walk
            continue;
        }
        result += walkResult;
    }
    return result / static_cast<double>(wpp_);
}

template<typename ScalarType, int DIM>
void
WoStSolver<ScalarType, DIM>::solve(const std::vector<Vector<DIM>>& points, std::vector<ScalarType>& results)
{
    const int numPoints = static_cast<int>(points.size());
    results.resize(numPoints);

    std::atomic<int> completedCount{ 0 };
    int lastReportedPct = -1;

    spdlog::info("Solving {} points at {} spp with {} threads", numPoints, wpp_, omp_get_max_threads());

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto threadRng = Solver<ScalarType, DIM>::makeRngFromSeed(this->seed_, 1 + tid);
        std::uniform_real_distribution<double> threadDist01(0.0, 1.0);

#pragma omp for schedule(dynamic)
        for (int i = 0; i < numPoints; ++i) {
            ScalarType sum(0.0);
            fcpw::Interaction<DIM> interaction;
            for (int w = 0; w < wpp_; ++w) {
                auto walkResult =
                  walkWoSt(*this->scene_, points[i], maxSteps_, epsilon_, rMin_, threadRng, threadDist01, interaction);
                if (walkResult.isNaN()) {
                    --w; // retry this walk
                    continue;
                }
                sum += walkResult;
            }
            results[i] = sum / static_cast<double>(wpp_);

            int done = completedCount.fetch_add(1, std::memory_order_relaxed) + 1;

            if (tid == 0) {
                int pct = done * 100 / numPoints;
                int pctStep = pct / 10;
                int lastPctStep = lastReportedPct / 10;
                if (pctStep != lastPctStep) {
                    spdlog::info("Progress: {}%", pct);
                    lastReportedPct = pct;
                }
            }
        }

        if (tid == 0)
            spdlog::info("Progress: 100%");
    }
}

template class WOS_SOLVER_API WoStSolver<Scalar<1>, 2>;
template class WOS_SOLVER_API WoStSolver<Scalar<1>, 3>;
template class WOS_SOLVER_API WoStSolver<Scalar<3>, 2>;
template class WOS_SOLVER_API WoStSolver<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
