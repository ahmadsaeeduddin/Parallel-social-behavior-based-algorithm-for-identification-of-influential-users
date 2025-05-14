
# PSAIIM Parallel Implementation

This repository contains multiple implementations of the **Parallel Social Behavior-Based Algorithm for Identification of Influential Users in Social Network (PSAIIM)**. The algorithm is based on the research paper:

**"Parallel social behavior-based algorithm for identification of influential users in social network"**  
by Wassim Mnasri et al.

The algorithm has been implemented using different parallel programming approaches to analyze and optimize its performance over large-scale social network graphs.

---

## 📁 Project Structure

Each subfolder in this project corresponds to a different implementation of the algorithm:

---

### 🔹 `Sequential/`

Contains the baseline serial (non-parallelized) implementation of PSAIIM.  
This version is useful for comparing performance with parallel implementations.

**Compilation:**
```bash
g++ -O2 -o out serial_main.cpp
```

---

### 🔹 `MPI/`

Implements the distributed memory version using **MPI**.  
This version divides the graph and spreads the computation across multiple MPI processes (ideal for clusters or multi-node systems).

**Compilation:**
```bash
mpic++ -o out main.cpp
```

**Execution Example:**
```bash
mpirun -np 4 --hostfile hosts.txt ./out
```

---

### 🔹 `OpenMP+MPI/`

This hybrid implementation combines **MPI** and **OpenMP**:
- MPI is used for inter-process communication
- OpenMP is used for intra-process multithreading

This version is optimized for running on multi-core multi-node systems.

**Compilation:**
```bash
mpic++ -fopenmp -o out main.cpp
```

**Execution Example:**
```bash
export OMP_NUM_THREADS=2
mpirun -np 4 --hostfile hosts.txt ./out
```

---

### 🔹 `OpenMP+MPI+METIS/`

The most advanced version of the implementation.  
This version uses:
- **METIS** for pre-partitioning the input graph
- **MPI** for process-level distribution
- **OpenMP** for shared-memory parallelism within each process

It is ideal for large, irregular graphs and demonstrates best scalability.

**Compilation:**
```bash
mpic++ -fopenmp -o out main.cpp -lmetis
```

---

## 📂 Dataset

All implementations are designed to work with the **Higgs Twitter social network dataset**.  
You should create a folder named:

```
Dataset-Higgs-Twitter/
```

inside each of the following directories:

- `MPI/`
- `OpenMP+MPI/`
- `OpenMP+MPI+METIS/`
- `Sequential/`

Each folder must contain the dataset files specific to that version, such as:
- `.edgelist` files
- METIS partition files (for the METIS version)
- ID mapping files (if needed)

---

### 🔗 Dataset Source

You can download the **Higgs Twitter social network dataset** from Stanford SNAP:

> 📥 [https://snap.stanford.edu/data/higgs-twitter.html](https://snap.stanford.edu/data/higgs-twitter.html)

Place the contents inside a folder named `Dataset-Higgs-Twitter` in each version folder as explained above.

---

## ⚙ Performance Profiling (VTune)

To analyze performance using **Intel VTune Profiler**, the project supports:

- Hotspot analysis
- Threading analysis
- MPI profiling (single-rank or per-rank)

VTune output folders (like `vtune_serial/`, `vtune_mpi/`, etc.) are **not committed** to the repo.

---

## 📄 .gitignore Guidelines

This project excludes the following via `.gitignore`:

```gitignore
# Ignore dataset in all folders
**/Dataset-Higgs-Twitter/

# Ignore VTune result directories
vtune_*/
**/vtune_*/

# Ignore binary files and core dumps
*.out
core
```

---

## 👤 Author

This codebase was developed as part of a semester project for the course **Parallel and Distributed Computing**.

The project explores performance tradeoffs of different parallel paradigms and showcases how domain knowledge (graph structure) and tools (e.g., METIS, VTune) can be combined to optimize complex real-world algorithms.
