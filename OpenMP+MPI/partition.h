#pragma once

#include <vector>
#include "node.h"    // your Node definition with out[], index, lowlink, compID, level, type, onStack

/// Compute SCC + CAC partitioning on the FOLLOW-graph stored in `graph`.
/// After this call:
///   - graph[v].compID is a 0-based component ID
///   - graph[v].type   is SCC or CAC
///   - graph[v].level  is its depth in the component DAG
#include <stack>
#include <queue>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <omp.h> // Include OpenMP header

void computeSCC_CAC(std::vector<Node>& graph) {
    const int N = (int)graph.size();

    // 1) Build FOLLOW adjacency and reverse adjacency
    std::vector<std::vector<int>> adj(N), radj(N);
    #pragma omp parallel for schedule(dynamic)
    for (int u = 0; u < N; ++u) {
        for (auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v >= 0 && v < N) {
                #pragma omp critical
                {
                    adj[u].push_back(v);
                    radj[v].push_back(u);
                }
            }
        }
    }

    // 2) First pass: DFS to compute finish order
    std::vector<char> seen(N, 0);
    std::vector<int> order;
    order.reserve(N);
    std::stack<int> dfs;

    #pragma omp parallel
    {
        std::vector<int> local_order;
        local_order.reserve(N);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            if (!seen[i]) {
                std::stack<int> local_dfs;
                local_dfs.push(i << 1); // even=enter, odd=exit
                while (!local_dfs.empty()) {
                    int x = local_dfs.top();
                    local_dfs.pop();
                    int v = x >> 1;
                    bool exiting = x & 1;
                    if (exiting) {
                        local_order.push_back(v);
                    } else if (!seen[v]) {
                        seen[v] = 1;
                        local_dfs.push((v << 1) | 1);
                        for (int w : adj[v]) {
                            if (!seen[w])
                                local_dfs.push(w << 1);
                        }
                    }
                }
            }
        }

        #pragma omp critical
        order.insert(order.end(), local_order.begin(), local_order.end());
    }

    // 3) Second pass: process in reverse finish order on radj
    std::vector<int> compID(N, -1);
    int C = 0;
    for (int idx = N - 1; idx >= 0; --idx) {
        int v = order[idx];
        if (compID[v] >= 0) continue;
        // collect one SCC
        dfs.push(v);
        compID[v] = C;
        std::vector<int> members;
        while (!dfs.empty()) {
            int u = dfs.top();
            dfs.pop();
            members.push_back(u);
            for (int w : radj[u]) {
                if (compID[w] < 0) {
                    compID[w] = C;
                    dfs.push(w);
                }
            }
        }
        C++;
    }

    // Initialize Node fields
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        graph[v].compID = compID[v];
        graph[v].type = Node::SCC; // mark all multi-node or singleton for now
        graph[v].level = 0;
    }

    // 4) Identify CACs: any component with size==1 and no self-loop
    std::vector<int> compSize(C, 0);
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        #pragma omp atomic
        compSize[compID[v]]++;
    }

    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        int c = compID[v];
        if (compSize[c] == 1) {
            graph[v].type = Node::CAC;
        }
    }

    // 5) Build component DAG edges (no duplicates)
    std::vector<std::vector<int>> cadj(C);
    #pragma omp parallel for
    for (int u = 0; u < N; ++u) {
        int cu = compID[u];
        for (int v : adj[u]) {
            int cv = compID[v];
            if (cu != cv) {
                #pragma omp critical
                cadj[cu].push_back(cv);
            }
        }
    }

    #pragma omp parallel for
    for (auto& nbrs : cadj) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    // 6) Compute topological levels via BFS from all sources
    std::vector<int> indeg(C, 0);
    for (int u = 0; u < C; ++u)
        for (int v : cadj[u]) indeg[v]++;

    std::queue<int> q;
    for (int u = 0; u < C; ++u)
        if (indeg[u] == 0) q.push(u);

    std::vector<int> clevel(C, 0);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v : cadj[u]) {
            clevel[v] = std::max(clevel[v], clevel[u] + 1);
            if (--indeg[v] == 0) q.push(v);
        }
    }

    // 7) Assign levels back to nodes
    #pragma omp parallel for
    for (int v = 0; v < N; ++v) {
        graph[v].level = clevel[compID[v]];
    }

    // 8) CAC‐merge: for each singleton comp, try merge one‐level‐down CAC neighbors
    std::vector<std::vector<int>> membersOf(C);
    for (int v = 0; v < N; ++v)
        membersOf[compID[v]].push_back(v);

    for (int c = 0; c < C; ++c) {
        if (membersOf[c].size() != 1) continue; // only singletons
        int u = membersOf[c][0];
        if (graph[u].type != Node::CAC) continue;
        int lvl = graph[u].level;
        // find CAC neighbors one level below
        for (auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v < 0 || v >= N) continue;
            int cv = compID[v];
            if (graph[v].type == Node::CAC && clevel[cv] == lvl - 1) {
                // merge component cv into c: redirect compID
                for (int w : membersOf[cv]) {
                    compID[w] = c;
                    graph[w].compID = c;
                }
                // no need to merge edges—levels remain consistent
                break;
            }
        }
    }

    // Done
    std::cout << "[LOG] computeSCC_CAC: " << C << " initial SCCs, "
              << "CAC-merged singleton chains\n";
}
