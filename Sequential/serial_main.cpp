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
using namespace std;

int main() {
    auto total_start = chrono::high_resolution_clock::now();

    constexpr uint32_t MAX_TEST_NODES = 10000;
    constexpr size_t N = MAX_TEST_NODES;

    cout << "[LOG] Initializing graph with " << N << " nodes\n";
    vector<Node> graph(N);

    cout << "[LOG] Generating interest vectors\n";
    generateInterestVectors(graph);
    cout << "[LOG] Interest vector generation complete.\n";

    cout << "[LOG] Loading graph layers\n";
    loadGraph(graph);
    cout << "[LOG] Graph load complete. Ready for SCC/CAC partitioning.\n";

    cout << "[LOG] Computing SCC/CAC partitioning...\n";
    computeSCC_CAC(graph);
    cout << "[LOG] Partitioning complete.\n";

    int maxCID = 0;
    for (auto& n : graph)
        if (n.compID > maxCID) maxCID = n.compID;
    cout << "[LOG] Found " << (maxCID + 1) << " components\n";

    cout << "[LOG] Computing influence power...\n";
    array<double, NUM_LAYERS> alpha = {0.0, 0.50, 0.15, 0.35}; // RETWEET, MENTION, REPLY weights
    double damping = 0.85;
    vector<double> IP = computeInfluencePower(graph, alpha, damping);
    cout << "[LOG] Influence power computation complete.\n";

    cout << "[LOG] Selecting seed candidates (Algorithm 6)...\n";
    vector<int> I_star = selectSeedCandidates(graph, IP);
    cout << "[LOG] Seed candidate selection complete.\n";

    cout << "[LOG] Final influential seed selection (Algorithm 7)...\n";
    vector<int> INF = selectFinalSeeds(graph, IP, I_star);

    vector<pair<int, double>> topSeeds;
    for (int u : INF) {
        topSeeds.emplace_back(u, IP[u]);
    }
    sort(topSeeds.begin(), topSeeds.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    cout << "[LOG] Top 10 final seed nodes (by IP):\n" << INF.size() <<endl<<endl<<endl;
    for (int i = 0; i < min(10, (int)topSeeds.size()); ++i) {
        cout << "  Node " << topSeeds[i].first << " → IP = " << setprecision(6) << topSeeds[i].second << "\n";
    }

    auto total_end = chrono::high_resolution_clock::now();
    cout << "[LOG] All stages complete.\n";
    cout << "[LOG] Total execution time: " << chrono::duration<double>(total_end - total_start).count() << " seconds.\n";
    return 0;
}
