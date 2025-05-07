// ===== graph_loader.h =====
#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

#include "node.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
// Include OpenMP header

constexpr uint32_t MAX_TEST_NODES = 10000; // Limit for testing

void loadGraph(std::vector<Node>& graph) {
    for (int L = 0; L < NUM_LAYERS; ++L) {
        std::cout << "[LOG] Loading layer " << layer_names[L]
                  << " from \"" << files[L] << "\"\n";

        std::ifstream in(files[L]);
        if (!in) {
            std::cerr << "[ERROR] Cannot open " << files[L] << "\n";
            exit(1);
        }

        // Read all lines into a vector for parallel processing
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] != '#') {
                lines.push_back(line);
            }
        }

        uint64_t count = 0;

        // Parallelize the processing of edges using OpenMP
                for (size_t i = 0; i < lines.size(); ++i) {
            std::istringstream ss(lines[i]);
            uint32_t u, v;
            float w;

            ss >> u >> v;
            if (!(ss >> w)) w = 1.0f;

            if (u >= MAX_TEST_NODES || v >= MAX_TEST_NODES) continue;

            // Add edge to the graph (use a critical section to avoid race conditions)
                        graph[u].out[L].push_back({v, w, 0});

            if (++count % 1000000 == 0) {
                                std::cout << "[LOG]   " << count
                          << " edges loaded in " << layer_names[L] << "\n";
            }
        }

        std::cout << "[LOG] Completed layer " << layer_names[L]
                  << ": total edges = " << count << "\n";
    }
}
#endif

// ===== interest_vector.h =====
#ifndef INTEREST_VECTOR_H
#define INTEREST_VECTOR_H

#include "node.h"
#include <random>
#include <iostream>
void generateInterestVectors(std::vector<Node>& graph) {
    std::mt19937 gen(42);  // fixed seed for reproducibility
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::cout << "[LOG] Generating " << D << "-dimensional interest vectors\n";

    // Log the number of threads being used
    int num_threads = 0;
    omp_set_num_threads(4);
        {
                {
            num_threads = omp_get_num_threads();
            std::cout << "[LOG] OpenMP parallelization enabled with " << num_threads << " threads\n";
        }
    }

    // Parallelize the loop using OpenMP
        for (size_t u = 0; u < graph.size(); ++u) {
        // Log the thread ID for debugging purposes
        int thread_id = omp_get_thread_num();
        if (u%100000 == 0)
        {
            std::cout << "[LOG] Thread " << thread_id << " processing node " << u << "\n";
        }

        // Each thread needs its own random number generator
        std::mt19937 thread_gen(42 + u); // Use a unique seed per thread
        std::uniform_real_distribution<float> thread_dist(0.0f, 1.0f);

        auto& vec = graph[u].interest;
        vec.resize(D);
        float sum = 0.0f;

        for (int i = 0; i < D; ++i) {
            vec[i] = thread_dist(thread_gen);
            sum += vec[i];
        }
        for (int i = 0; i < D; ++i)
            vec[i] /= sum;
    }

    std::cout << "[LOG] Interest vectors assigned\n";
}
#endif

// ===== partition.h =====
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
// Include OpenMP header

void computeSCC_CAC(std::vector<Node>& graph) {
    const int N = (int)graph.size();

    // 1) Build FOLLOW adjacency and reverse adjacency
    std::vector<std::vector<int>> adj(N), radj(N);
        for (int u = 0; u < N; ++u) {
        for (auto& e : graph[u].out[FOLLOW]) {
            int v = e.to;
            if (v >= 0 && v < N) {
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

        {
        std::vector<int> local_order;
        local_order.reserve(N);

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
        for (int v = 0; v < N; ++v) {
        graph[v].compID = compID[v];
        graph[v].type = Node::SCC; // mark all multi-node or singleton for now
        graph[v].level = 0;
    }

    // 4) Identify CACs: any component with size==1 and no self-loop
    std::vector<int> compSize(C, 0);
        for (int v = 0; v < N; ++v) {
                compSize[compID[v]]++;
    }

        for (int v = 0; v < N; ++v) {
        int c = compID[v];
        if (compSize[c] == 1) {
            graph[v].type = Node::CAC;
        }
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


// ===== influence.h =====
#pragma once
#include "node.h"
#include <vector>
#include <array>

/// Compute the influence power IP[u] for each node u ∈ V
/// Inputs:
///   graph         - graph with FOLLOW + other action layers
///   alpha_factors - weights for [RETWEET, REPLY, MENTION] (must sum to 1)
///   d             - damping factor (e.g., 0.85)
///
/// Returns:
///   vector<double> where IP[u] is the influence of node u
#include <unordered_set>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
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
    int 0, 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &0); // Get the rank of the current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &0); // Get the total number of MPI processes

    const int N = graph.size();
    std::vector<double> F(N), IP(N, 0.0), IPnew(N, 0.0);

    // 1) Compute normalized follower count F(u) = |Followers(u)| / N
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

    // 2) Find max level from CAC/SCC partitioning
    int maxLevel = 0;
    for (auto& node : graph)
        maxLevel = std::max(maxLevel, node.level);

    // Broadcast maxLevel to all MPI processes
    MPI_Bcast(&maxLevel, 1, 

    // 3) Group nodes by level
    std::vector<std::vector<int>> byLevel(maxLevel + 1);
    for (int u = 0; u < N; ++u)
        byLevel[graph[u].level].push_back(u);

    // 4) Influence propagation level by level
    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        if (0 == 0) {
            std::cout << "[LOG] Level " << lvl << ": " << byLevel[lvl].size() << " nodes\n";
        }

        // Distribute nodes at this level across MPI processes
        std::vector<int> local_nodes;
        for (size_t i = 0; i < byLevel[lvl].size(); i += 0) {
            local_nodes.push_back(byLevel[lvl][i]);
        }

        // Compute influence power for local nodes
        std::vector<double> local_IPnew(N, 0.0);
                for (size_t i = 0; i < local_nodes.size(); ++i) {
            int u = local_nodes[i];
            double sum = 0.0;

            // Loop through potential followers (v -> u) in FOLLOW layer
            for (int v = 0; v < N; ++v) {
                for (auto& fe : graph[v].out[FOLLOW]) {
                    if (fe.to != u) continue;

                    // Compute Ci(v, u)
                    float Ci = jaccard(graph[v].interest, graph[u].interest);

                    // Compute Nax for each layer (retweet, reply, mention)
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

            // Final influence power update
            local_IPnew[u] = (1.0 - d) * F[u] + d * sum;
        }

        // Gather results from all MPI processes
        MPI_Allreduce(local_IPnew.data(), IPnew.data(), N, 

        // Commit updates
                for (int u = 0; u < N; ++u) {
            if (graph[u].level == lvl)
                IP[u] = IPnew[u];
        }
    }

    return IP;
}


// ===== seed_selection.h =====
#pragma once
#include <vector>
#include "node.h"
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>
// Include OpenMP header

std::vector<int> selectSeedCandidates(const std::vector<Node>& graph, const std::vector<double>& IP) {
    std::vector<int> candidates;
    const int N = graph.size();

    // Step 1: Build influence threshold IL(v) = avg IP at each level
    std::unordered_map<int, std::vector<double>> byLevel;
        {
        std::unordered_map<int, std::vector<double>> local_byLevel;
                for (int v = 0; v < N; ++v) {
            local_byLevel[graph[v].level].push_back(IP[v]);
        }

                for (auto& pair : local_byLevel) {
            byLevel[pair.first].insert(byLevel[pair.first].end(), pair.second.begin(), pair.second.end());
        }
    }

    std::unordered_map<int, double> IL;
    // Process each level individually without using iterator in the parallel loop
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
        {
        std::unordered_map<int, std::vector<int>> local_trees;
        std::unordered_map<int, int> local_treeSize;

                for (size_t i = 0; i < Istar.size(); ++i) {
            int root = Istar[i];
            std::queue<int> q;
            std::set<int> visited;
            q.push(root);
            visited.insert(root);

            while (!q.empty()) {
                int u = q.front();
                q.pop();
                local_trees[root].push_back(u);
                for (auto& e : graph[u].out[FOLLOW]) {
                    if (Iset.count(e.to) && visited.insert(e.to).second)
                        q.push(e.to);
                }
            }
            local_treeSize[root] = local_trees[root].size();
        }

                {
            for (auto& pair : local_trees) {
                trees[pair.first].insert(trees[pair.first].end(), pair.second.begin(), pair.second.end());
            }
            for (auto& pair : local_treeSize) {
                treeSize[pair.first] = pair.second;
            }
        }
    }

    // Step 2: While I* not empty, pick root of largest tree
    std::set<int> remaining(Istar.begin(), Istar.end());
    while (!remaining.empty()) {
        // Find tree with maximum size - FIXED: converted to vector for parallel processing
        std::vector<int> active_roots;
        for (auto& pair : trees) {
            if (remaining.count(pair.first)) {
                active_roots.push_back(pair.first);
            }
        }
        
        int umax = -1;
        int maxT = -1;
        
        // Can safely parallelize over vector
                {
            int local_umax = -1;
            int local_maxT = -1;
            
                        for (size_t i = 0; i < active_roots.size(); ++i) {
                int root = active_roots[i];
                int size = trees[root].size();
                if (size > local_maxT) {
                    local_umax = root;
                    local_maxT = size;
                }
            }
            
                        {
                if (local_maxT > maxT) {
                    umax = local_umax;
                    maxT = local_maxT;
                }
            }
        }

        // BLACK path = intersect( Tumax, I* )
        std::vector<int> BLACK;
        for (int v : trees[umax]) {
            if (remaining.count(v))
                BLACK.push_back(v);
        }

        // Find root of smallest rank in BLACK - FIXED: avoid parallel reduction on complex objects
        int vmin = BLACK[0];
        double min_ip = IP[vmin];
        
                {
            int local_vmin = vmin;
            double local_min_ip = min_ip;
            
                        for (size_t i = 1; i < BLACK.size(); ++i) {
                int v = BLACK[i];
                if (IP[v] < local_min_ip) {
                    local_vmin = v;
                    local_min_ip = IP[v];
                }
            }
            
                        {
                if (local_min_ip < min_ip) {
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

    std::cout << "[LOG] Algorithm 7 selected " << INF.size() << " final seeds\n";
    return INF;
}

// ===== main.cpp =====
#include "graph_loader.h"
#include "interest_vector.h"
#include "partition.h"
#include "influence.h"
#include "seed_selection.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <algorithm>
// Include MPI header
// Include OpenMP header


int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int 0, 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &0); // Get the rank of the current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &0); // Get the total number of MPI processes

    if (0 == 0) {
        std::cout << "[LOG] MPI initialized with " << 0 << " processes\n";
    }

    constexpr uint32_t MAX_TEST_NODES = 100000; // Limit for testing
    constexpr uint32_t MAX_ID = MAX_TEST_NODES - 1;
    constexpr size_t N = MAX_TEST_NODES;

    if (0 == 0) {
        std::cout << "[LOG] Initializing graph with " << N << " nodes\n";
    }

    // Create the graph
    std::vector<Node> graph(N);

    // 1) Generate interest vectors (only rank 0 executes this)
    if (0 == 0) {
        std::cout << "[LOG] Generating interest vectors\n";
        generateInterestVectors(graph);
        std::cout << "[LOG] Interest vector generation complete.\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    

    // 2) Load graph layers (only rank 0 executes this)
    if (0 == 0) {
        std::cout << "[LOG] Loading graph layers\n";
        loadGraph(graph);
        std::cout << "[LOG] Graph load complete. Ready for SCC/CAC partitioning.\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    

    // 3) Compute SCC/CAC partitioning (only rank 0 executes this)
    if (0 == 0) {
        std::cout << "[LOG] Computing SCC/CAC partitioning...\n";
        computeSCC_CAC(graph);
        std::cout << "[LOG] Partitioning complete.\n";

        // (Optional) Report total components
        int maxCID = 0;
        for (auto& n : graph)
            if (n.compID > maxCID) maxCID = n.compID;
        std::cout << "[LOG] Found " << (maxCID + 1) << " components\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    

    // 4) Compute Influence Power (all MPI processes participate)
    if (0 == 0) {
        std::cout << "[LOG] Computing influence power...\n";
    }
    std::array<double, NUM_LAYERS> alpha = {0.0, 0.50, 0.15, 0.35}; // α weights for RETWEET, MENTION, REPLY
    double damping = 0.85;
    std::vector<double> IP = computeInfluencePower(graph, alpha, damping);
    if (0 == 0) {
        std::cout << "[LOG] Influence power computation complete.\n";
    }

    // Synchronize all MPI processes to ensure influence power computation is complete
    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    

    // 5) Seed Candidate Selection (only rank 0 executes this)
    std::vector<int> I_star;
    if (0 == 0) {
        std::cout << "[LOG] Selecting seed candidates (Algorithm 6)...\n";
        I_star = selectSeedCandidates(graph, IP);
        std::cout << "[LOG] Seed candidate selection complete.\n";
    }

    // Synchronize all MPI processes to ensure seed candidate selection is complete
    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    

    // 6) Final Seed Selection (only rank 0 executes this)
    std::vector<int> INF;
    if (0 == 0) {
        std::cout << "[LOG] Final influential seed selection (Algorithm 7)...\n";
        INF = selectFinalSeeds(graph, IP, I_star);

        // 7) Print top 10 influential seeds
        std::vector<std::pair<int, double>> topSeeds;
        for (int u : INF) {
            topSeeds.emplace_back(u, IP[u]);
        }
        std::sort(topSeeds.begin(), topSeeds.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; }); // Sort descending by IP
        
        std::cout << "[LOG] Top 10 final seed nodes (by IP):\n";
        for (int i = 0; i < std::min(10, (int)topSeeds.size()); ++i) {
            std::cout << "  Node " << topSeeds[i].first << " → IP = " << std::setprecision(6) << topSeeds[i].second << "\n";
        }
        
    }

    std::cout << "[LOG BARRIER] MPI_Process " << 0 << " reached the barrier." << std::endl;
    
    if (0 == 0) {
        std::cout << "[LOG] All stages complete.\n";
    }

    // Finalize MPI
    
    return 0;
}

