#pragma once
#include <vector>
#include "node.h"
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>

std::vector<int> selectSeedCandidates(const std::vector<Node>& graph, const std::vector<double>& IP) {
    std::vector<int> candidates;
    const int N = graph.size();

    // Step 1: Build influence threshold IL(v) = avg IP at each level
    std::unordered_map<int, std::vector<double>> byLevel;
    for (int v = 0; v < N; ++v) {
        byLevel[graph[v].level].push_back(IP[v]);
    }

    std::unordered_map<int, double> IL;
    for (auto& pair : byLevel) {
        int level = pair.first;
        std::vector<double>& values = pair.second;

        double sum = 0;
        for (size_t i = 0; i < values.size(); ++i) {
            sum += values[i];
        }
        IL[level] = sum / std::max(1.0, (double)values.size());
    }

    // Step 2: Iterate each node to check candidate criteria
    for (int v = 0; v < N; ++v) {
        int L = graph[v].level;
        double IPL = IL[L];
        double IPLp1 = IL.count(L + 1) ? IL[L + 1] : 0.0;

        if ((IPL - IPLp1) > IPLp1 && IP[v] > IPL) {
            candidates.push_back(v);
        }
    }

    std::cout << "[LOG] Algorithm 6 selected " << candidates.size() << " seed candidates\n";
    return candidates;
}

std::vector<int> selectFinalSeeds(const std::vector<Node>& graph,
                                  const std::vector<double>& IP,
                                  const std::vector<int>& Istar)
{
    std::vector<int> INF;
    std::set<int> Iset(Istar.begin(), Istar.end());
    std::unordered_map<int, std::vector<int>> trees;
    std::unordered_map<int, int> treeSize;

    // Step 1: Build Influence-BFS Trees
    for (size_t i = 0; i < Istar.size(); ++i) {
        int root = Istar[i];
        std::queue<int> q;
        std::set<int> visited;
        q.push(root);
        visited.insert(root);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            trees[root].push_back(u);
            for (auto& e : graph[u].out[FOLLOW]) {
                if (Iset.count(e.to) && visited.insert(e.to).second)
                    q.push(e.to);
            }
        }
        treeSize[root] = trees[root].size();
    }

    // Step 2: While I* not empty, pick root of largest tree
    std::set<int> remaining(Istar.begin(), Istar.end());
    while (!remaining.empty()) {
        std::vector<int> active_roots;
        for (auto& pair : trees) {
            if (remaining.count(pair.first)) {
                active_roots.push_back(pair.first);
            }
        }

        int umax = -1;
        int maxT = -1;
        for (size_t i = 0; i < active_roots.size(); ++i) {
            int root = active_roots[i];
            int size = trees[root].size();
            if (size > maxT) {
                umax = root;
                maxT = size;
            }
        }

        std::vector<int> BLACK;
        for (int v : trees[umax]) {
            if (remaining.count(v))
                BLACK.push_back(v);
        }

        int vmin = BLACK[0];
        double min_ip = IP[vmin];
        for (size_t i = 1; i < BLACK.size(); ++i) {
            int v = BLACK[i];
            if (IP[v] < min_ip) {
                vmin = v;
                min_ip = IP[v];
            }
        }

        INF.push_back(vmin);

        for (int v : BLACK)
            remaining.erase(v);
        remaining.erase(vmin);
    }

    std::cout << "[LOG] Algorithm 7 selected " << INF.size() << " final seeds\n";
    return INF;
}
