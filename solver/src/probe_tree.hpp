#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <unordered_set>
#include <vector>

#include "solver/probe.hpp"

#include "core/vector.hpp"

WOS_NAMESPACE_OPEN_SCOPE

namespace detail {

template<typename ScalarType, int DIM>
struct ProbeTree
{
    using ProbeType = Probe<ScalarType, DIM>;

    struct LeafEntry
    {
        int probeIndex;
        double ratioMin; // min ||p-c||/R_c for any p in this leaf
        double ratioMax; // max ||p-c||/R_c for any p in this leaf
    };

    struct Node
    {
        std::array<std::unique_ptr<Node>, 1 << DIM> children;
        std::vector<LeafEntry>
          entries; // In leaf entries, we only store probes that have the potential to be an optimal probe
        std::vector<int> containingProbes; // Here we store all probes that contain any point in this node, for fast
                                           // lookup of containing probes
        Vector<DIM> center;
        double halfSize = 0.0;
        bool isLeaf = true;
        int depth = 0;
    };

    std::unique_ptr<Node> root;
    const std::vector<ProbeType>* probes_ = nullptr;
    double alphaRec_;
    double alphaWalk_;

    ProbeTree(double alphaRec, double alphaWalk)
      : alphaRec_(alphaRec)
      , alphaWalk_(alphaWalk)
    {
    }

    static constexpr int maxDepth = 8;
    static constexpr int maxPerLeaf = 4;

    // ==== Queries ====

    const ProbeType* findOptimal(const std::vector<ProbeType>& probes, const Vector<DIM>& p, double alpha) const
    {
        if (!root)
            return nullptr;

        const Node* node = root.get();

        while (!node->isLeaf) {
            int ci = 0;
            for (int d = 0; d < DIM; ++d)
                if (p[d] >= node->center[d])
                    ci |= (1 << d);
            node = node->children[ci].get();
            if (!node)
                return nullptr;
        }

        const ProbeType* bestProbe = nullptr;
        double bestNormDist = std::numeric_limits<double>::max();

        for (const auto& e : node->entries) {
            if (e.ratioMin >= alpha)
                continue;

            const auto& probe = probes[e.probeIndex];
            double d = (p - probe.center).norm();
            double nd = d / probe.radius;
            if (nd < alpha && nd < bestNormDist) {
                bestNormDist = nd;
                bestProbe = &probe;
            }
        }

        return bestProbe;
    }

    const std::vector<const ProbeType*> findContaining(const std::vector<ProbeType>& probes,
                                                       const Vector<DIM>& p,
                                                       double alpha) const
    {
        if (!root)
            return {};

        const Node* node = root.get();

        while (!node->isLeaf) {
            int ci = 0;
            for (int d = 0; d < DIM; ++d)
                if (p[d] >= node->center[d])
                    ci |= (1 << d);
            node = node->children[ci].get();
            if (!node)
                return {};
        }

        std::vector<const ProbeType*> result;

        for (int probeIdx : node->containingProbes) {
            const auto& probe = probes[probeIdx];
            if (probe.contains(p, alpha))
                result.push_back(&probe);
        }

        return result;
    }

    // ==== Build ====

    void build(const std::vector<ProbeType>& probes, const Vector<DIM>& bboxMin, const Vector<DIM>& bboxMax)
    {
        probes_ = &probes;

        if (probes.empty())
            return;

        root = std::make_unique<Node>();
        root->center = (bboxMin + bboxMax) * 0.5;
        root->halfSize = (bboxMax - bboxMin).maxCoeff() * 0.5;
        root->isLeaf = true;
        root->depth = 0;

        for (int i = 0; i < static_cast<int>(probes.size()); ++i)
            insert(i, root.get());
    }

  private:
    // ==== Helpers ====

    static Vector<DIM> closestPoint(const Node& node, const Vector<DIM>& p)
    {
        Vector<DIM> cp;
        for (int d = 0; d < DIM; ++d)
            cp[d] = std::clamp(p[d], node.center[d] - node.halfSize, node.center[d] + node.halfSize);
        return cp;
    }

    static Vector<DIM> farthestPoint(const Node& node, const Vector<DIM>& p)
    {
        Vector<DIM> fp;
        for (int d = 0; d < DIM; ++d) {
            double lo = node.center[d] - node.halfSize;
            double hi = node.center[d] + node.halfSize;
            fp[d] = (p[d] < (lo + hi) * 0.5) ? hi : lo;
        }
        return fp;
    }

    // ==== Insert ====

    void insert(int probeIndex, Node* node, bool containmentOnly = false)
    {
        const auto& probe = (*probes_)[probeIndex];
        double effRadius = std::max(alphaRec_, alphaWalk_) * probe.radius;

        Vector<DIM> cp = closestPoint(*node, probe.center);
        double nearestDist = (cp - probe.center).norm();
        if (nearestDist > effRadius)
            return;

        if (node->isLeaf)
            insertToLeaf(probeIndex, probe, node, nearestDist, containmentOnly);
        else
            for (auto& child : node->children)
                if (child)
                    insert(probeIndex, child.get(), containmentOnly);
    }

    void insertToLeaf(int probeIndex, const ProbeType& probe, Node* node, double nearestDist, bool containmentOnly)
    {
        node->containingProbes.push_back(probeIndex);
        if (containmentOnly)
            return;

        Vector<DIM> fp = farthestPoint(*node, probe.center);
        double farthestDist = (fp - probe.center).norm();

        double newRatioMin = nearestDist / probe.radius;
        double newRatioMax = farthestDist / probe.radius;

        bool dominated = false;
        for (auto it = node->entries.begin(); it != node->entries.end();) {
            if (newRatioMin >= it->ratioMax) {
                dominated = true;
                break;
            } else if (newRatioMax <= it->ratioMin) {
                it = node->entries.erase(it);
            } else {
                ++it;
            }
        }

        if (!dominated) {
            node->entries.push_back({ probeIndex, newRatioMin, newRatioMax });

            if (static_cast<int>(node->entries.size()) > maxPerLeaf && node->depth < maxDepth)
                subdivide(node);
        }
    }

    // ==== Subdivide ====

    void subdivide(Node* node)
    {
        node->isLeaf = false;
        double childHalf = node->halfSize * 0.5;

        for (int i = 0; i < (1 << DIM); ++i) {
            auto child = std::make_unique<Node>();
            child->halfSize = childHalf;
            child->depth = node->depth + 1;
            child->isLeaf = true;
            for (int d = 0; d < DIM; ++d) {
                double sign = ((i >> d) & 1) ? 1.0 : -1.0;
                child->center[d] = node->center[d] + sign * childHalf;
            }
            node->children[i] = std::move(child);
        }

        std::vector<LeafEntry> oldEntries;
        oldEntries.swap(node->entries);
        std::vector<int> oldContaining;
        oldContaining.swap(node->containingProbes);

        for (const auto& e : oldEntries)
            insert(e.probeIndex, node);

        std::unordered_set<int> oldEntriesSet;
        for (const auto& e : oldEntries)
            oldEntriesSet.insert(e.probeIndex);
        for (int probeIdx : oldContaining)
            if (oldEntriesSet.find(probeIdx) == oldEntriesSet.end())
                insert(probeIdx, node, true);
    }
};

} // namespace detail

WOS_NAMESPACE_CLOSE_SCOPE
