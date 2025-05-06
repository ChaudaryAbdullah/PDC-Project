Parallel SSSP Update Framework for Dynamic Networks

A scalable, platform-independent parallel implementation of a Single-Source Shortest Path (SSSP) update algorithm for large-scale dynamic networks, leveraging MPI, OpenMP, OpenCL, and METIS. This project is developed as part of the Parallel and Distributed Computing (PDC) course at FAST-NUCES, Islamabad.
📖 Project Overview
Real-world networks, such as social graphs, communication systems, and biological networks, are large (millions of vertices, billions of edges) and dynamic, with frequent edge insertions and deletions. Traditional SSSP algorithms like Dijkstra’s are inefficient for dynamic graphs, requiring costly recomputation. Inspired by Khanda et al. (2022), this project implements a two-phase parallel algorithm template that incrementally updates SSSP trees, avoiding full recomputation. Our implementation scales across distributed-memory clusters and heterogeneous systems (CPUs and GPUs) using:

METIS: For graph partitioning to minimize inter-node communication.
MPI: For inter-node parallelism and coordination.
OpenMP: For intra-node shared-memory parallelism.
OpenCL: For GPU-accelerated edge relaxation.

The project includes serial and parallel SSSP implementations, an OpenCL kernel for GPU offloading, and a Python script for performance visualization.

👥 Team
Abdullah Shakir
Messam Raza 
Arban Arfan 

🚀 Features
Dynamic SSSP Updates: Efficiently processes edge insertions and deletions without recomputing the entire SSSP tree.
Two-Phase Algorithm:
Step 1: Identifies affected subgraphs using parallel edge processing.
Step 2: Iteratively updates affected vertices with lock-free relaxations.


Hybrid Parallelization:
METIS for load-balanced graph partitioning.
MPI for distributed-memory coordination.
OpenMP for multi-core CPU parallelism.
OpenCL for GPU acceleration.


Performance Visualization: Python script (plotGraph.py) generates bar charts comparing execution times across different MPI process counts.
Scalability: Supports large graphs (up to 16M vertices, 250M edges) and dynamic update batches.
Robustness: Handles malformed inputs, negative weights, and self-loops with detailed error reporting.

🛠️ Implementation Details
The repository contains the following key components:
C++ Implementations

Serial SSSP (serial_execution.cpp):
Implements Dijkstra’s algorithm for static and dynamic SSSP.
Processes graph and update files, applies edge changes, and saves results.
Command: g++ -std=c++11 serial_execution.cpp -o serial_sssp && ./serial_sssp sample_graph.txt sample_updates.txt 10000 output.txt


Parallel SSSP (main.cpp, graph.cpp, sssp.cpp, utils.cpp, opencl_utils.cpp):
Distributed implementation using MPI for inter-node parallelism.
OpenMP for intra-node parallelism with dynamic scheduling.
METIS for graph partitioning to minimize communication overhead.
OpenCL for GPU-accelerated edge relaxation.
Command: mpicxx -O3 -march=native -funroll-loops -fopenmp -DCL_TARGET_OPENCL_VERSION=200 -o sssp main.cpp graph.cpp utils.cpp sssp.cpp opencl_utils.cpp -I. -L/usr/local/lib -lOpenCL -lmetis && mpirun --use-hwthread-cpus --bind-to core:overload-allowed -np 4 ./sssp sample_graph.txt sample_updates.txt 10000 output.txt --openmp --opencl


OpenCL Kernel

Edge Relaxation (relax_edges.cl):
Optimized kernel for parallel edge relaxation on GPUs.
Uses atomic operations for thread-safe distance updates.
Integrated via opencl_utils.cpp for seamless CPU-GPU coordination.


Python Visualization
Performance Plotting (plotGraph.py):
Executes the parallel SSSP program with varying MPI processes (2, 3, 4).
Generates a bar chart (sssp_performance.png) of execution times.
Command: python plotGraph.py


📊 Algorithm Workflow
Based on Khanda et al. (2022), the algorithm operates in two phases:
Identify Affected Subgraph:
Processes edge insertions/deletions in parallel (MPI/OpenMP/OpenCL).
Marks affected vertices and subtrees using boolean flags (Affected, Affected_Del).
Deletions disconnect subtrees (set distances to ∞); insertions tentatively relax paths.

Update Affected Subgraph:
Iteratively relaxes marked vertices’ edges until convergence.
Uses lock-free, asynchronous relaxations to minimize synchronization overhead.
MPI exchanges boundary flags to propagate cross-partition effects.

Key Data Structures:

Adjacency list for graph representation.
SSSP tree with parent pointers, distances, and affected flags.
Arrays for inserted (Ins_k) and deleted (Del_k) edges.

📈 Performance Highlights
Speedup: Achieves up to 8.5× speedup over Gunrock (GPU) and 5× over Galois (CPU) for moderate update sizes (<75% deletions).
Scalability: Tested on graphs with ~16M vertices and 250M edges.
Load Balancing: METIS partitioning and OpenMP dynamic scheduling ensure even workload distribution.
Heterogeneous Computing: Seamlessly integrates CPU (OpenMP) and GPU (OpenCL) processing.

🧰 Prerequisites
OS: Linux (Ubuntu recommended) or Windows with WSL2 for MPI/OpenCL.
Compilers:
g++ (GCC 9.x+) with C++17 support.
mpicxx (OpenMPI 4.x or MPICH 3.x) for parallel compilation.


Libraries:
MPI: OpenMPI or MPICH for distributed parallelism.
OpenMP: For shared-memory parallelism (included with GCC).
OpenCL: Version 2.0+ for GPU acceleration (e.g., NVIDIA, AMD, or Intel drivers).
METIS: Version 5.1.0+ for graph partitioning (libmetis).


Python: 3.6+ with matplotlib for visualization (pip install matplotlib).
Hardware:
Multi-core CPU (recommended: 8+ cores).
GPU with OpenCL support (optional but recommended).
Sufficient RAM for large graphs (16GB+ recommended).

