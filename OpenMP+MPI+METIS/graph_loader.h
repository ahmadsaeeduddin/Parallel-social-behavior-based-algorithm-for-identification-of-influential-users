#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

#include "node.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <omp.h> // Include OpenMP header

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
        #pragma omp parallel for schedule(dynamic) reduction(+:count)
        for (size_t i = 0; i < lines.size(); ++i) {
            std::istringstream ss(lines[i]);
            uint32_t u, v;
            float w;

            ss >> u >> v;
            if (!(ss >> w)) w = 1.0f;

            if (u >= MAX_TEST_NODES || v >= MAX_TEST_NODES) continue;

            // Add edge to the graph (use a critical section to avoid race conditions)
            #pragma omp critical
            graph[u].out[L].push_back({v, w, 0});

            if (++count % 1000000 == 0) {
                #pragma omp critical
                std::cout << "[LOG]   " << count
                          << " edges loaded in " << layer_names[L] << "\n";
            }
        }

        std::cout << "[LOG] Completed layer " << layer_names[L]
                  << ": total edges = " << count << "\n";
    }
}
#endif