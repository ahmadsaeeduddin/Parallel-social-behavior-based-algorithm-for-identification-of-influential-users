#include "serial_graph_loader.h"
#include "serial_interest_vector.h"
#include "serial_partition.h"
#include "serial_influence.h"
#include "serial_seed_selection.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <algorithm>
#include <chrono>

int main() {
    auto total_start = std::chrono::high_resolution_clock::now();

    constexpr uint32_t MAX_TEST_NODES = 10000;
    constexpr size_t N = MAX_TEST_NODES;

    std::cout << "[LOG] Initializing graph with " << N << " nodes\n";
    std::vector<Node> graph(N);

    std::cout << "[LOG] Generating interest vectors\n";
    generateInterestVectors(graph);
    std::cout << "[LOG] Interest vector generation complete.\n";

    std::cout << "[LOG] Loading graph layers\n";
    loadGraph(graph);
    std::cout << "[LOG] Graph load complete. Ready for SCC/CAC partitioning.\n";

    std::cout << "[LOG] Computing SCC/CAC partitioning...\n";
    computeSCC_CAC(graph);
    std::cout << "[LOG] Partitioning complete.\n";

    int maxCID = 0;
    for (auto& n : graph)
        if (n.compID > maxCID) maxCID = n.compID;
    std::cout << "[LOG] Found " << (maxCID + 1) << " components\n";

    std::cout << "[LOG] Computing influence power...\n";
    std::array<double, NUM_LAYERS> alpha = {0.0, 0.50, 0.15, 0.35}; // RETWEET, MENTION, REPLY weights
    double damping = 0.85;
    std::vector<double> IP = computeInfluencePower(graph, alpha, damping);
    std::cout << "[LOG] Influence power computation complete.\n";

    std::cout << "[LOG] Selecting seed candidates (Algorithm 6)...\n";
    std::vector<int> I_star = selectSeedCandidates(graph, IP);
    std::cout << "[LOG] Seed candidate selection complete.\n";

    std::cout << "[LOG] Final influential seed selection (Algorithm 7)...\n";
    std::vector<int> INF = selectFinalSeeds(graph, IP, I_star);

    std::vector<std::pair<int, double>> topSeeds;
    for (int u : INF) {
        topSeeds.emplace_back(u, IP[u]);
    }
    std::sort(topSeeds.begin(), topSeeds.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    std::cout << "[LOG] Top 10 final seed nodes (by IP):\n";
    for (int i = 0; i < std::min(10, (int)topSeeds.size()); ++i) {
        std::cout << "  Node " << topSeeds[i].first << " → IP = " << std::setprecision(6) << topSeeds[i].second << "\n";
    }

    auto total_end = std::chrono::high_resolution_clock::now();
    std::cout << "[LOG] All stages complete.\n";
    std::cout << "[LOG] Total execution time: " << std::chrono::duration<double>(total_end - total_start).count() << " seconds.\n";
    return 0;
}
