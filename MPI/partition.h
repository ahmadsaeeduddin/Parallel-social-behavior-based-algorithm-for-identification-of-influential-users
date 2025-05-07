#pragma once

#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <cassert>
#include <iostream>
#include "node.h"

void computeSCC_CAC(std::vector<Node>& graph) {
    const int N = (int)graph.size();

    // 1) Build FOLLOW adjacency and reverse adjacency
    std::vector<std::vector<int>> adj(N), radj(N);
    for (int u = 0; u < N; ++u) {
        for (auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v >= 0 && v < N) {
                adj[u].push_back(v);
                radj[v].push_back(u);
            }
        }
    }

    // 2) First pass: DFS to compute finish order
    std::vector<char> seen(N, 0);
    std::vector<int> order;
    order.reserve(N);
    std::stack<int> dfs;

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
                    order.push_back(v);
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

    // 3) Second pass: process in reverse finish order on radj
    std::vector<int> compID(N, -1);
    int C = 0;
    for (int idx = N - 1; idx >= 0; --idx) {
        int v = order[idx];
        if (compID[v] >= 0) continue;
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
    for (int v = 0; v < N; ++v) {
        graph[v].compID = compID[v];
        graph[v].type = Node::SCC;
        graph[v].level = 0;
    }

    // 4) Identify CACs: components with size 1
    std::vector<int> compSize(C, 0);
    for (int v = 0; v < N; ++v)
        compSize[compID[v]]++;

    for (int v = 0; v < N; ++v) {
        int c = compID[v];
        if (compSize[c] == 1)
            graph[v].type = Node::CAC;
    }

    // 5) Build component DAG edges (no duplicates)
    std::vector<std::vector<int>> cadj(C);
    for (int u = 0; u < N; ++u) {
        int cu = compID[u];
        for (int v : adj[u]) {
            int cv = compID[v];
            if (cu != cv) {
                cadj[cu].push_back(cv);
            }
        }
    }

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
    for (int v = 0; v < N; ++v)
        graph[v].level = clevel[compID[v]];

    // 8) CAC‐merge: for each singleton comp, try merge one‐level‐down CAC neighbors
    std::vector<std::vector<int>> membersOf(C);
    for (int v = 0; v < N; ++v)
        membersOf[compID[v]].push_back(v);

    for (int c = 0; c < C; ++c) {
        if (membersOf[c].size() != 1) continue;
        int u = membersOf[c][0];
        if (graph[u].type != Node::CAC) continue;
        int lvl = graph[u].level;
        for (auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v < 0 || v >= N) continue;
            int cv = compID[v];
            if (graph[v].type == Node::CAC && clevel[cv] == lvl - 1) {
                for (int w : membersOf[cv]) {
                    compID[w] = c;
                    graph[w].compID = c;
                }
                break;
            }
        }
    }

    std::cout << "[LOG] computeSCC_CAC: " << C << " initial SCCs, "
              << "CAC-merged singleton chains\n";
}
