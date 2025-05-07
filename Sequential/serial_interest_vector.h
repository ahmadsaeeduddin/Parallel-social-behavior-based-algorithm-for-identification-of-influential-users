#ifndef INTEREST_VECTOR_H
#define INTEREST_VECTOR_H

#include "node.h"
#include <random>
#include <iostream>
using namespace std;

void generateInterestVectors(vector<Node>& graph) {
    cout << "[LOG] Generating interest vectors\n";

    for (int u = 0; u < graph.size(); ++u) {
        vector<float> interests(D);
        float sum = 0.0f;

        // Generate D random values
        for (int i = 0; i < D; ++i) {
            interests[i] = static_cast<float>(rand()) / RAND_MAX;
            sum += interests[i];
        }

        // Normalize
        for (int i = 0; i < D; ++i) {
            interests[i] /= sum;
        }

        graph[u].interest = interests;
    }

    cout << "[LOG] Interest vectors assigned\n";
}


#endif
