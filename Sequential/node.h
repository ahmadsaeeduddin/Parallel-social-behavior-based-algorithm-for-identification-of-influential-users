#ifndef NODE_H
#define NODE_H

#include <vector>
#include <array>
#include <cstdint>
#include <string>
<<<<<<< HEAD
using namespace std;

enum Layer : uint8_t { FOLLOW = 0, RETWEET = 1, MENTION = 2, REPLY = 3, NUM_LAYERS };

static const char* layer_names[NUM_LAYERS] = { "FOLLOW", "RETWEET", "MENTION", "REPLY" };

static const array<string, NUM_LAYERS> files = {
=======

enum Layer : uint8_t { FOLLOW = 0, RETWEET = 1, MENTION = 2, REPLY = 3, NUM_LAYERS };
static const char* layer_names[NUM_LAYERS] = { "FOLLOW", "RETWEET", "MENTION", "REPLY" };
static const std::array<std::string, NUM_LAYERS> files = {
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d
    "Dataset-Higgs-Twitter/higgs-social_network.edgelist",
    "Dataset-Higgs-Twitter/higgs-retweet_network.edgelist",
    "Dataset-Higgs-Twitter/higgs-mention_network.edgelist",
    "Dataset-Higgs-Twitter/higgs-reply_network.edgelist"
};

static constexpr int D = 10;

struct Edge {
    uint32_t to;
    float    weight;
    uint32_t timestamp;
};

struct Node {
    // adjacency lists for each layer
<<<<<<< HEAD
    array<vector<Edge>, NUM_LAYERS> out;

    // interest vector
    vector<float> interest;
=======
    std::array<std::vector<Edge>, NUM_LAYERS> out;

    // interest vector
    std::vector<float> interest;
>>>>>>> 39fc27665bc0d2d95f0a79a90378e87e4dbe239d

    // --- SCC/CAC fields ---
    int index   = -1;             // Tarjan discovery index
    int lowlink = -1;             // Tarjan lowlink
    bool onStack = false;         // is it on Tarjan’s stack?
    int compID  = -1;             // final component ID
    int level   = 0;              // level in the component DAG
    enum Type { UNCLASSIFIED, SCC, CAC } type = UNCLASSIFIED;
};

#endif // NODE_H
