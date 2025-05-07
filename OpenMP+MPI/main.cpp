#include "graph_loader.h"
#include "interest_vector.h"
#include "partition.h"
#include "influence.h"
#include "seed_selection.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <algorithm>
#include <mpi.h> // Include MPI header
#include <omp.h> // Include OpenMP header
using namespace std;

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    double start_time = 0.0, end_time = 0.0;
    start_time = MPI_Wtime();

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank); // Get the rank of the current MPI process
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size); // Get the total number of MPI processes

    if (mpi_rank == 0)
    {
        cout << "[LOG] MPI initialized with " << mpi_size << " processes\n";
    }

    const uint32_t MAX_TEST_NODES = 10000; // Limit for testing
    const uint32_t MAX_ID = MAX_TEST_NODES - 1;
    const size_t N = MAX_TEST_NODES;

    if (mpi_rank == 0)
    {
        cout << "[LOG] Initializing graph with " << N << " nodes\n";
    }

    // Create the graph
    vector<Node> graph(N);

    // 1) Generate interest vectors (only rank 0 executes this)
    if (mpi_rank == 0)
    {
        cout << "[LOG] Generating interest vectors\n";
        generateInterestVectors(graph);
        cout << "[LOG] Interest vector generation complete.\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 2) Load graph layers (only rank 0 executes this)
    if (mpi_rank == 0)
    {
        cout << "[LOG] Loading graph layers\n";
        loadGraph(graph);
        cout << "[LOG] Graph load complete. Ready for SCC/CAC partitioning.\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    // cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 3) Compute SCC/CAC partitioning (only rank 0 executes this)
    if (mpi_rank == 0)
    {
        cout << "[LOG] Computing SCC/CAC partitioning...\n";
        computeSCC_CAC(graph);
        cout << "[LOG] Partitioning complete.\n";

        // (Optional) Report total components
        int maxCID = 0;
        for (auto &n : graph)
            if (n.compID > maxCID)
                maxCID = n.compID;
        cout << "[LOG] Found " << (maxCID + 1) << " components\n";
    }

    // Synchronize all MPI processes to ensure rank 0 completes before others proceed
    cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 4) Compute Influence Power (all MPI processes participate)
    if (mpi_rank == 0)
    {
        cout << "[LOG] Computing influence power...\n";
    }
    array<double, NUM_LAYERS> alpha = {0.0, 0.50, 0.15, 0.35}; // α weights for RETWEET, MENTION, REPLY
    double damping = 0.85;
    vector<double> IP = computeInfluencePower(graph, alpha, damping);
    if (mpi_rank == 0)
    {
        cout << "[LOG] Influence power computation complete.\n";
    }

    // Synchronize all MPI processes to ensure influence power computation is complete
    cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 5) Seed Candidate Selection (only rank 0 executes this)
    vector<int> I_star;
    if (mpi_rank == 0)
    {
        cout << "[LOG] Selecting seed candidates (Algorithm 6)...\n";
        I_star = selectSeedCandidates(graph, IP);
        cout << "[LOG] Seed candidate selection complete.\n";
    }

    // Synchronize all MPI processes to ensure seed candidate selection is complete
    cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    // 6) Final Seed Selection (only rank 0 executes this)
    vector<int> INF;
    if (mpi_rank == 0)
    {
        cout << "[LOG] Final influential seed selection (Algorithm 7)...\n";
        INF = selectFinalSeeds(graph, IP, I_star);

        // 7) Print top 10 influential seeds
        vector<pair<int, double>> topSeeds;
        for (int u : INF)
        {
            topSeeds.emplace_back(u, IP[u]);
        }
        sort(topSeeds.begin(), topSeeds.end(),
             [](const auto &a, const auto &b)
             { return a.second > b.second; }); // Sort descending by IP

        cout << "[LOG] Top 10 final seed nodes (by IP):\n";
        for (int i = 0; i < min(10, (int)topSeeds.size()); ++i)
        {
            cout << "  Node " << topSeeds[i].first << " → IP = " << setprecision(6) << topSeeds[i].second << "\n";
        }
    }

    cout << "[LOG BARRIER] MPI_Process " << mpi_rank << " reached the barrier." << endl;
    MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes reach this point
    end_time = MPI_Wtime();

    if (mpi_rank == 0)
    {
        cout << "[LOG] All stages complete.\n";
        cout << "[LOG] Total execution time: " << (end_time - start_time) << " seconds.\n";
    }

    MPI_Finalize();
    return 0;
}