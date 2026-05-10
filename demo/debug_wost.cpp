#include "core/scalar.hpp"
#include "core/vector.hpp"
#include "scene/obj_scene.hpp"
#include "solver/wost_solver.hpp"
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using Scalar3d = WOS::Scalar<3>;
using Vec2 = WOS::Vector<2>;

int
main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: debug_wost <scene_dir> <x> <y> [num_walks] [seed] [max_steps]\n";
        return 1;
    }

    std::string scenePath = argv[1];
    double px = std::stod(argv[2]);
    double py = std::stod(argv[3]);
    int numWalks = (argc >= 5) ? std::stoi(argv[4]) : 1000;
    uint64_t seed = (argc >= 6) ? std::stoull(argv[5]) : 42;
    int maxSteps = (argc >= 7) ? std::stoi(argv[6]) : 1024;

    spdlog::set_level(spdlog::level::info);
    spdlog::info("Loading scene: {}", scenePath);

    WOS::ObjScene<Scalar3d, 2> scene(scenePath);
    Vec2 p(px, py);

    spdlog::info("Test point: ({}, {})", px, py);
    spdlog::info("ROI check: {}", scene.isInsideROI(p));

    // Check boundary distances
    {
        fcpw::Interaction<2> inter;
        if (scene.findClosestDirichletBoundary(p, inter))
            spdlog::info("Closest Dirichlet boundary: dist={}, pos=[{},{}]",
                         inter.d,
                         inter.p[0],
                         inter.p[1]);
        else
            spdlog::info("No Dirichlet boundary found!");

        if (scene.findClosestNeumannBoundary(p, inter))
            spdlog::info("Closest Neumann boundary: dist={}, pos=[{},{}]",
                         inter.d,
                         inter.p[0],
                         inter.p[1]);
        else
            spdlog::info("No Neumann boundary found!");
    }

    // Check source at point
    spdlog::info("Source at point: [{}, {}, {}]",
                 scene.source(p)[0],
                 scene.source(p)[1],
                 scene.source(p)[2]);

    spdlog::info("Running {} WoSt walks with seed {}, maxSteps={}...",
                 numWalks,
                 seed,
                 maxSteps);

    WOS::WoStSolver<Scalar3d, 2> solver(scene, seed);
    solver.setMaxSteps(maxSteps);
    solver.setWalksPerPixel(1);

    std::vector<Scalar3d> results;
    std::vector<double> rVals, gVals, bVals;
    results.reserve(numWalks);
    rVals.reserve(numWalks);
    gVals.reserve(numWalks);
    bVals.reserve(numWalks);

    int negCount = 0;
    int zeroCount = 0;
    double rMin = std::numeric_limits<double>::max();
    double rMax = std::numeric_limits<double>::lowest();

    for (int i = 0; i < numWalks; ++i) {
        auto val = solver.solve(p);
        results.push_back(val);
        rVals.push_back(val[0]);
        gVals.push_back(val[1]);
        bVals.push_back(val[2]);

        if (val[0] < 0.0 || val[1] < 0.0 || val[2] < 0.0)
            negCount++;
        if (val[0] == 0.0 && val[1] == 0.0 && val[2] == 0.0)
            zeroCount++;

        rMin = std::min(rMin, val[0]);
        rMax = std::max(rMax, val[0]);

        if ((i + 1) % 500 == 0)
            spdlog::info("  {}/{} walks done, neg={}, zero={}, min={}, max={}",
                         i + 1,
                         numWalks,
                         negCount,
                         zeroCount,
                         rMin,
                         rMax);
    }

    // Statistics
    auto stats = [](const std::vector<double>& v) {
        double sum = 0.0, sumSq = 0.0;
        double mn = std::numeric_limits<double>::max();
        double mx = std::numeric_limits<double>::lowest();
        for (auto x : v) {
            sum += x;
            sumSq += x * x;
            mn = std::min(mn, x);
            mx = std::max(mx, x);
        }
        double mean = sum / v.size();
        double var = sumSq / v.size() - mean * mean;
        return std::make_tuple(mn, mx, mean, std::sqrt(std::max(0.0, var)));
    };

    auto [r_min, r_max, r_mean, r_std] = stats(rVals);
    auto [g_min, g_max, g_mean, g_std] = stats(gVals);
    auto [b_min, b_max, b_mean, b_std] = stats(bVals);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n===== Results =====\n";
    std::cout << "Channel R: min=" << r_min << " max=" << r_max << " mean=" << r_mean
              << " std=" << r_std << "\n";
    std::cout << "Channel G: min=" << g_min << " max=" << g_max << " mean=" << g_mean
              << " std=" << g_std << "\n";
    std::cout << "Channel B: min=" << b_min << " max=" << b_max << " mean=" << b_mean
              << " std=" << b_std << "\n";
    std::cout << "Negative walks: " << negCount << " / " << numWalks << " ("
              << (100.0 * negCount / numWalks) << "%)\n";
    std::cout << "Zero walks: " << zeroCount << " / " << numWalks << " ("
              << (100.0 * zeroCount / numWalks) << "%)\n";

    // Show first 10 negative results
    std::cout << "\nFirst 10 negative sample values:\n";
    int shown = 0;
    for (int i = 0; i < numWalks && shown < 10; ++i) {
        if (results[i][0] < 0.0) {
            std::cout << "  walk " << i << ": [" << results[i][0] << ", " << results[i][1]
                      << ", " << results[i][2] << "]\n";
            shown++;
        }
    }

    return 0;
}
