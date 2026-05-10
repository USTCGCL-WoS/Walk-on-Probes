#pragma once

#include "core/math_defs.hpp"
#include "core/scalar.hpp"
#include "core/vector.hpp"

#include <vector>

WOS_NAMESPACE_OPEN_SCOPE

static constexpr int MaxWalkSteps = 2048;

template<typename ScalarType, int DIM>
struct WalkPath
{
    std::vector<Vector<DIM>> positions;
    std::vector<ScalarType> sourceContribs;
    std::vector<ScalarType> neumannContribs;
    std::vector<double> pdfs;                 // sampling PDF for each step's angle
    std::vector<int> probeIndices;            // -1 = fallback WoSt step
    ScalarType dirichlet = ScalarType::NaN(); // NaN if not hit Dirichlet boundary
    Vector<DIM> dirichletPos;
    int dirichletProbeIdx = -1;
    int count = 0; // number of recorded steps

    WalkPath()
    {
        positions.reserve(64);
        sourceContribs.reserve(64);
        neumannContribs.reserve(64);
        pdfs.reserve(64);
        probeIndices.reserve(64);
    }

    void clear()
    {
        positions.clear();
        sourceContribs.clear();
        neumannContribs.clear();
        pdfs.clear();
        probeIndices.clear();
        dirichlet = ScalarType::NaN();
        dirichletProbeIdx = -1;
        count = 0;
    }

    void recordStep(const Vector<DIM>& pos, double pdf, ScalarType S, ScalarType N, int probeIdx)
    {
        positions.push_back(pos);
        pdfs.push_back(pdf);
        sourceContribs.push_back(S);
        neumannContribs.push_back(N);
        probeIndices.push_back(probeIdx);
        ++count;
    }

    // O(M) reverse pass: estimates[0..count-1] for positions[0..count-1],
    // estimates[count] for dirichletPos (= D)
    void computeEstimates(std::vector<ScalarType>& estimates, int& estCount) const
    {
        estimates.resize(count + 1);
        ScalarType accum = dirichlet;
        estimates[count] = accum;
        for (int j = count - 1; j >= 0; --j) {
            accum += sourceContribs[j] + neumannContribs[j];
            estimates[j] = accum;
        }
        estCount = count + 1;
    }

    ScalarType startingEstimate() const
    {
        ScalarType r = dirichlet;
        for (int j = count - 1; j >= 0; --j)
            r += sourceContribs[j] + neumannContribs[j];
        return r;
    }
};

WOS_NAMESPACE_CLOSE_SCOPE
