#pragma once
#include <vector>
#include "node.h"
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
using namespace std;

// Include OpenMP header

vector<int> selectSeedCandidates(const vector<Node> &graph, const vector<double> &IP)
{
    vector<int> candidates;
    const int N = graph.size();

    // Step 1: Build influence threshold IL(v) = avg IP at each level
    unordered_map<int, vector<double>> byLevel;
    {
        unordered_map<int, vector<double>> local_byLevel;
        for (int v = 0; v < N; ++v)
        {
            local_byLevel[graph[v].level].push_back(IP[v]);
        }

        for (auto &pair : local_byLevel)
        {
            byLevel[pair.first].insert(byLevel[pair.first].end(), pair.second.begin(), pair.second.end());
        }
    }

    unordered_map<int, double> IL;
    // Process each level individually without using iterator in the parallel loop
    for (auto &pair : byLevel)
    {
        int level = pair.first;
        vector<double> &values = pair.second;

        double sum = 0;
        for (size_t i = 0; i < values.size(); ++i)
        {
            sum += values[i];
        }
        IL[level] = sum / max(1.0, (double)values.size());
    }

    // Step 2: Iterate each node to check candidate criteria
    for (int v = 0; v < N; ++v)
    {
        int L = graph[v].level;
        double IPL = IL[L];
        double IPLp1 = IL.count(L + 1) ? IL[L + 1] : 0.0;

        if ((IPL - IPLp1) > IPLp1 && IP[v] > IPL)
        {
            candidates.push_back(v);
        }
    }

    cout << "[LOG] Algorithm 6 selected " << candidates.size() << " seed candidates\n";
    return candidates;
}

vector<int> selectFinalSeeds(const vector<Node> &graph,
                             const vector<double> &IP,
                             const vector<int> &Istar)
{
    vector<int> INF;
    set<int> Iset(Istar.begin(), Istar.end());
    unordered_map<int, vector<int>> trees;
    unordered_map<int, int> treeSize;

    // Step 1: Build Influence-BFS Trees
    {
        unordered_map<int, vector<int>> local_trees;
        unordered_map<int, int> local_treeSize;

        for (size_t i = 0; i < Istar.size(); ++i)
        {
            int root = Istar[i];
            queue<int> q;
            set<int> visited;
            q.push(root);
            visited.insert(root);

            while (!q.empty())
            {
                int u = q.front();
                q.pop();
                local_trees[root].push_back(u);
                for (auto &e : graph[u].out[FOLLOW])
                {
                    if (Iset.count(e.to) && visited.insert(e.to).second)
                        q.push(e.to);
                }
            }
            local_treeSize[root] = local_trees[root].size();
        }

        {
            for (auto &pair : local_trees)
            {
                trees[pair.first].insert(trees[pair.first].end(), pair.second.begin(), pair.second.end());
            }
            for (auto &pair : local_treeSize)
            {
                treeSize[pair.first] = pair.second;
            }
        }
    }

    // Step 2: While I* not empty, pick root of largest tree
    set<int> remaining(Istar.begin(), Istar.end());
    while (!remaining.empty())
    {
        // Find tree with maximum size - FIXED: converted to vector for parallel processing
        vector<int> active_roots;
        for (auto &pair : trees)
        {
            if (remaining.count(pair.first))
            {
                active_roots.push_back(pair.first);
            }
        }

        int umax = -1;
        int maxT = -1;

        // Can safely parallelize over vector
        {
            int local_umax = -1;
            int local_maxT = -1;

            for (size_t i = 0; i < active_roots.size(); ++i)
            {
                int root = active_roots[i];
                int size = trees[root].size();
                if (size > local_maxT)
                {
                    local_umax = root;
                    local_maxT = size;
                }
            }

            {
                if (local_maxT > maxT)
                {
                    umax = local_umax;
                    maxT = local_maxT;
                }
            }
        }

        // BLACK path = intersect( Tumax, I* )
        vector<int> BLACK;
        for (int v : trees[umax])
        {
            if (remaining.count(v))
                BLACK.push_back(v);
        }

        // Find root of smallest rank in BLACK - FIXED: avoid parallel reduction on complex objects
        int vmin = BLACK[0];
        double min_ip = IP[vmin];

        {
            int local_vmin = vmin;
            double local_min_ip = min_ip;

            for (size_t i = 1; i < BLACK.size(); ++i)
            {
                int v = BLACK[i];
                if (IP[v] < local_min_ip)
                {
                    local_vmin = v;
                    local_min_ip = IP[v];
                }
            }

            {
                if (local_min_ip < min_ip)
                {
                    vmin = local_vmin;
                    min_ip = local_min_ip;
                }
            }
        }

        INF.push_back(vmin);

        // Remove BLACK + root from I*
        for (int v : BLACK)
            remaining.erase(v);
        remaining.erase(vmin);
    }

    cout << "[LOG] Algorithm 7 selected " << INF.size() << " final seeds\n";
    return INF;
}