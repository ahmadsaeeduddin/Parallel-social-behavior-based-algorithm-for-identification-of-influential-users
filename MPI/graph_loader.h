#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

#include "node.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

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

        std::string line;
        uint64_t count = 0;

        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;

            std::istringstream ss(line);
            uint32_t u, v;
            float w;

            ss >> u >> v;
            if (!(ss >> w)) w = 1.0f;

            if (u >= MAX_TEST_NODES || v >= MAX_TEST_NODES) continue;

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
