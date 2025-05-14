#ifndef GRAPH_LOADER_H
#define GRAPH_LOADER_H

#include "node.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
<<<<<<< HEAD
using namespace std;

const int MAX_TEST_NODES = 10000;

void loadGraph(vector<Node>& graph) {
    for (int L = 0; L < NUM_LAYERS; ++L) {
        cout << "[LOG] Loading layer " << layer_names[L] << " from \"" << files[L] << "\"\n";

        ifstream in(files[L]);
        if (!in) {
            cerr << "[ERROR] Cannot open " << files[L] << "\n";
=======
// Include OpenMP header

constexpr uint32_t MAX_TEST_NODES = 10000; // Limit for testing

void loadGraph(std::vector<Node>& graph) {
    for (int L = 0; L < NUM_LAYERS; ++L) {
        std::cout << "[LOG] Loading layer " << layer_names[L]
                  << " from \"" << files[L] << "\"\n";

        std::ifstream in(files[L]);
        if (!in) {
            std::cerr << "[ERROR] Cannot open " << files[L] << "\n";
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d
            exit(1);
        }

        // Read all lines into a vector for parallel processing
<<<<<<< HEAD
        vector<string> lines;
        string line;
        while (getline(in, line)) {
=======
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d
            if (!line.empty() && line[0] != '#') {
                lines.push_back(line);
            }
        }

<<<<<<< HEAD
        int count = 0;

        for (int i = 0; i < lines.size(); ++i) {
            istringstream ss(lines[i]);
=======
        uint64_t count = 0;

        // Parallelize the processing of edges using OpenMP
                for (size_t i = 0; i < lines.size(); ++i) {
            std::istringstream ss(lines[i]);
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d
            uint32_t u, v;
            float w;

            ss >> u >> v;
            if (!(ss >> w)) w = 1.0f;

            if (u >= MAX_TEST_NODES || v >= MAX_TEST_NODES) continue;

<<<<<<< HEAD
            graph[u].out[L].push_back({v, w, 0});

            if (++count % 10000 == 0) {
                cout << "[LOG]   " << count << " edges loaded in " << layer_names[L] << "\n";
            }
        }

        cout << "[LOG] Completed layer " << layer_names[L] << ": total edges = " << count << "\n";
=======
            // Add edge to the graph (use a critical section to avoid race conditions)
                        graph[u].out[L].push_back({v, w, 0});

            if (++count % 1000000 == 0) {
                                std::cout << "[LOG]   " << count
                          << " edges loaded in " << layer_names[L] << "\n";
            }
        }

        std::cout << "[LOG] Completed layer " << layer_names[L]
                  << ": total edges = " << count << "\n";
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d
    }
}
#endif