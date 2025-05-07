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
#include <omp.h>
#include <mpi.h>
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

inline void printProgress(int current, int total, int mpi_rank) {
    int barWidth = 40;
    float progress = float(current) / total;

    std::cout << "\r[Rank " << mpi_rank << "] [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
        std::cout << (i < pos ? "=" : " ");
    std::cout << "] " << int(progress * 100.0) << "%";
    std::cout.flush();
}


std::vector<double> computeInfluencePower(const std::vector<Node>& localGraph,
    const std::array<double, NUM_LAYERS>& alpha,
    double d,
    int N)
{
    const int localN = localGraph.size(); // Local node count

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get the rank of the current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get the total number of MPI processes

    //const int N = localGraph.size();
    std::vector<double> F(N), IP(N, 0.0), IPnew(N, 0.0);

    // 1) Compute normalized follower count F(u) = |Followers(u)| / N
    std::vector<int> followers(N, 0);
    #pragma omp parallel for
    for (int u = 0; u < N; ++u) {
        for (auto& e : localGraph[u].out[FOLLOW]) {
            if (e.to < N) {
                #pragma omp atomic
                followers[e.to]++;
            }
        }
    }
    #pragma omp parallel for
    for (int u = 0; u < N; ++u)
        F[u] = (double)followers[u] / N;

    // 2) Find max level from CAC/SCC partitioning
    int maxLevel = 0;
    for (auto& node : localGraph)
        maxLevel = std::max(maxLevel, node.level);

    // Broadcast maxLevel to all MPI processes
    MPI_Bcast(&maxLevel, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 3) Group nodes by level
    std::vector<std::vector<int>> byLevel(maxLevel + 1);
    for (int u = 0; u < N; ++u)
        byLevel[localGraph[u].level].push_back(u);

    // 4) Influence propagation level by level
    for (int lvl = 0; lvl <= maxLevel; ++lvl) {
        if (mpi_rank == 0) {
            std::cout << "[LOG] Level " << lvl << ": " << byLevel[lvl].size() << " nodes\n";
        }

        // Distribute nodes at this level across MPI processes
        std::vector<int> local_nodes;
        for (size_t i = mpi_rank; i < byLevel[lvl].size(); i += mpi_size) {
            local_nodes.push_back(byLevel[lvl][i]);
        }

        // Compute influence power for local nodes
        std::vector<double> local_IPnew(N, 0.0);
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < local_nodes.size(); ++i) {
            if (omp_get_thread_num() == 0 && i % 500 == 0) {
                #pragma omp critical
                printProgress(i, local_nodes.size(), mpi_rank);
            }

            int u = local_nodes[i];
            double sum = 0.0;

            // Loop through potential followers (v -> u) in FOLLOW layer
            for (int v = 0; v < N; ++v) {
                for (auto& fe : localGraph[v].out[FOLLOW]) {
                    if (fe.to != u) continue;

                    float Ci = jaccard(localGraph[v].interest, localGraph[u].interest);
                    std::array<int, NUM_LAYERS> Nax = {0};
                    for (auto& re : localGraph[v].out[RETWEET]) if (re.to == u) Nax[RETWEET]++;
                    for (auto& ce : localGraph[v].out[REPLY])   if (ce.to == u) Nax[REPLY]++;
                    for (auto& me : localGraph[v].out[MENTION]) if (me.to == u) Nax[MENTION]++;

                    double psi = 0.0;
                    for (int l = 1; l < NUM_LAYERS; ++l)
                        psi += alpha[l] * Ci * Nax[l];

                    int deg_v = std::max(1, (int)localGraph[v].out[FOLLOW].size());
                    sum += psi * (IP[v] / deg_v);
                }
            }

            local_IPnew[u] = (1.0 - d) * F[u] + d * sum;
        }


        // Gather results from all MPI processes
        MPI_Allreduce(local_IPnew.data(), IPnew.data(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Commit updates
        #pragma omp parallel for
        for (int u = 0; u < N; ++u) {
            if (localGraph[u].level == lvl)
                IP[u] = IPnew[u];
        }
        if (mpi_rank == 0) std::cout << std::endl; // Newline after progress bar
    }

    return IP;
}
