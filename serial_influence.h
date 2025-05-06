#pragma once
#include "node.h"
#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <algorithm>

static float jaccard(const std::vector<float>& a, const std::vector<float>& b) {
    float dot = 0, sum_a = 0, sum_b = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot += std::min(a[i], b[i]);
        sum_a += std::max(a[i], b[i]);
    }
    return sum_a > 0 ? dot / sum_a : 0.0f;
}

std::vector<double> computeInfluencePower(const std::vector<Node>& graph, const std::array<double, NUM_LAYERS>& alpha, double d)
{
    const int N = graph.size();
    std::vector<double> F(N), IP(N, 0.0), IPnew(N, 0.0);

    // 1) Compute normalized follower count F(u)
    std::vector<int> followers(N, 0);
    for (int u = 0; u < N; ++u) {
        for (auto& e : graph[u].out[FOLLOW]) {
            if (e.to < N) {
                followers[e.to]++;
            }
        }
    }
    for (int u = 0; u < N; ++u)
        F[u] = (double)followers[u] / N;

    // 2) Find max level from SCC/CAC partitioning
    int maxLevel = 0;
    for (auto& node : graph)
        maxLevel = std::max(maxLevel, node.level);

    // 3) Group nodes by level
    std::vector<std::vector<int>> byLevel(maxLevel + 1);
    for (int u = 0; u < N; ++u)
        byLevel[graph[u].level].push_back(u);

    // 4) Influence propagation level by level
    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        std::cout << "[LOG] Level " << lvl << ": " << byLevel[lvl].size() << " nodes\n";

        for (int i = 0; i < byLevel[lvl].size(); ++i) {
            int u = byLevel[lvl][i];
            double sum = 0.0;

            // Loop through potential followers (v → u)
            for (int v = 0; v < N; ++v) {
                for (auto& fe : graph[v].out[FOLLOW]) {
                    if (fe.to != u) continue;

                    // Compute Ci(v, u)
                    float Ci = jaccard(graph[v].interest, graph[u].interest);

                    // Compute Nax
                    std::array<int, NUM_LAYERS> Nax = {0};
                    for (auto& re : graph[v].out[RETWEET])
                        if (re.to == u) Nax[RETWEET]++;
                    for (auto& ce : graph[v].out[REPLY])
                        if (ce.to == u) Nax[REPLY]++;
                    for (auto& me : graph[v].out[MENTION])
                        if (me.to == u) Nax[MENTION]++;

                    // Compute ψ(v,u)
                    double psi = 0.0;
                    for (int l = 1; l < NUM_LAYERS; ++l)
                        psi += alpha[l] * Ci * Nax[l];

                    int deg_v = std::max(1, (int)graph[v].out[FOLLOW].size());
                    sum += psi * (IP[v] / deg_v);
                }
            }

            // Update influence score
            IPnew[u] = (1.0 - d) * F[u] + d * sum;
        }

        // Commit updates for this level
        for (int u = 0; u < N; ++u) {
            if (graph[u].level == lvl)
                IP[u] = IPnew[u];
        }
    }

    return IP;
}
