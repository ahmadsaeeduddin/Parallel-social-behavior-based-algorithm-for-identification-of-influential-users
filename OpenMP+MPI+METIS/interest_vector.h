#ifndef INTEREST_VECTOR_H
#define INTEREST_VECTOR_H

#include "node.h"
#include <random>
#include <iostream>
#include <omp.h>

void generateInterestVectors(std::vector<Node>& graph) {
    std::mt19937 gen(42);  // fixed seed for reproducibility
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::cout << "[LOG] Generating " << D << "-dimensional interest vectors\n";

    // Log the number of threads being used
    int num_threads = 0;
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            std::cout << "[LOG] OpenMP parallelization enabled with " << num_threads << " threads\n";
        }
    }

    // Parallelize the loop using OpenMP
    #pragma omp parallel for
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