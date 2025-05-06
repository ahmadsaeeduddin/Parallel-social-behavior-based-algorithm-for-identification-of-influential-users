#ifndef INTEREST_VECTOR_H
#define INTEREST_VECTOR_H

#include "node.h"
#include <random>
#include <iostream>

void generateInterestVectors(std::vector<Node>& graph) {
    std::mt19937 gen(42);  // fixed seed for reproducibility
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::cout << "[LOG] Generating " << D << "-dimensional interest vectors\n";

    for (size_t u = 0; u < graph.size(); ++u) {
        if (u % 100000 == 0) {
            std::cout << "[LOG] Processing node " << u << "\n";
        }

        // Use a unique seed per node for reproducibility
        std::mt19937 thread_gen(42 + u);
        std::uniform_real_distribution<float> thread_dist(0.0f, 1.0f);

        auto& vec = graph[u].interest;
        vec.resize(D);
        float sum = 0.0f;

        for (int i = 0; i < D; ++i) {
            vec[i] = thread_dist(thread_gen);
            sum += vec[i];
        }

        for (int i = 0; i < D; ++i) {
            vec[i] /= sum;
        }
    }

    std::cout << "[LOG] Interest vectors assigned\n";
}

#endif
