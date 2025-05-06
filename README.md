⚡ Parallel SSSP Update Framework for Dynamic Networks

A scalable, platform-independent parallel implementation of the Single-Source Shortest Path (SSSP) update algorithm for large-scale dynamic networks, leveraging MPI, OpenMP, OpenCL, and METIS.
Developed for the Parallel and Distributed Computing (PDC) course at FAST-NUCES, Islamabad.

📖 Project Overview

Real-world networks — like social graphs, communication systems, and biological networks — are large and dynamic, with millions of vertices and frequent edge insertions/deletions. Traditional SSSP algorithms like Dijkstra’s are inefficient in such environments due to costly recomputation.

Inspired by Khanda et al. (2022), this project implements a two-phase parallel algorithm that incrementally updates the SSSP tree, avoiding full recomputation. Our implementation scales across distributed-memory clusters and heterogeneous systems using:

🧩 METIS for graph partitioning

🖧 MPI for inter-node parallelism

🧵 OpenMP for shared-memory parallelism

🖥️ OpenCL for GPU-accelerated edge relaxation

👥 Team

Abdullah Shakir

Messam Raza

Arban Arfan

🚀 Features
✅ Dynamic SSSP Updates
Efficiently handles edge insertions and deletions without full recomputation.

⚙️ Two-Phase Algorithm

Identify Affected Subgraph: Detects affected vertices via parallel edge processing.

Update Phase: Iteratively relaxes affected vertices using lock-free, asynchronous updates.

🧵 Hybrid Parallelization

METIS: Load-balanced partitioning

MPI: Distributed-memory coordination

OpenMP: CPU parallelism

OpenCL: GPU acceleration

📊 Performance Visualization

Python script generates bar charts comparing execution times across MPI process counts.

🏗️ Scalability

Supports graphs with up to 16M vertices and 250M edges.

🛡️ Robustness

Handles malformed inputs, negative weights, and self-loops with detailed error handling.

🛠️ Implementation Details

📁 Components

1. Serial SSSP (serial_execution.cpp)
   
Implements Dijkstra’s algorithm for static and dynamic graphs.

Applies edge updates and stores results.


g++ -std=c++11 serial_execution.cpp -o serial_sssp
./serial_sssp sample_graph.txt sample_updates.txt 10000 output.txt

2. Parallel SSSP (main.cpp, graph.cpp, sssp.cpp, utils.cpp, opencl_utils.cpp)
   
MPI + OpenMP + METIS + OpenCL integration for full hybrid parallelism.


mpicxx -O3 -march=native -funroll-loops -fopenmp -DCL_TARGET_OPENCL_VERSION=200 \
-o sssp main.cpp graph.cpp utils.cpp sssp.cpp opencl_utils.cpp -I. \
-L/usr/local/lib -lOpenCL -lmetis

mpirun --use-hwthread-cpus --bind-to core:overload-allowed -np 4 \
./sssp sample_graph.txt sample_updates.txt 10000 output.txt --openmp --opencl

3. OpenCL Kernel (relax_edges.cl)
   
Optimized GPU kernel using atomic operations for safe edge relaxations.

4. Python Visualization (plotGraph.py)

Benchmarks performance across 2–4 MPI processes.

Generates sssp_performance.png.

📊 Algorithm Workflow

🔁 Two-Phase Update (based on Khanda et al. 2022)

Phase 1: Identify Affected Subgraph

Processes edge insertions/deletions in parallel.

Flags affected vertices (Affected, Affected_Del).

Handles deletions by disconnecting subtrees; handles insertions by tentative relaxation.

Phase 2: Update Affected Subgraph

Lock-free, asynchronous edge relaxations until convergence.

MPI communicates boundary updates for cross-partition consistency.

📚 Key Data Structures

Adjacency List: Efficient graph representation.

SSSP Tree: Stores distances, parents, and update flags.

Edge Arrays: Ins_k, Del_k for dynamic updates.

⚙️ Performance Highlights

🚀 Speedup: Up to 8.5× over Gunrock (GPU) and 5× over Galois (CPU) for moderate updates.

🌐 Scalability: Tested on graphs with ~16M vertices, 250M edges.

⚖️ Load Balancing: Achieved using METIS and OpenMP’s dynamic scheduling.

🔀 Heterogeneous Computing: Unified CPU (OpenMP) + GPU (OpenCL) processing pipeline.

🧰 Prerequisites

🖥️ OS: 
Linux (Ubuntu recommended) or Windows with WSL2

🔧 Compilers

g++ (GCC 9.x+) with C++17 support

mpicxx (OpenMPI 4.x / MPICH 3.x)

📦 Required Libraries

MPI (OpenMPI or MPICH)

OpenMP (bundled with GCC)

OpenCL (v2.0+)

METIS (v5.1.0+)

🐍 Python

Python 3.6+

Install dependencies with:

pip install matplotlib

💻 Hardware

Multi-core CPU (8+ cores recommended)

OpenCL-compatible GPU (optional but recommended)

RAM: 16 GB+ for large graphs
