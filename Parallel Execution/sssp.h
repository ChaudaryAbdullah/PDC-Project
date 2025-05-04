
#ifndef SSSP_H
#define SSSP_H

#include "graph.h"
#include "opencl_utils.h"
#include <vector>
#include <mpi.h>

class SSSP
{
public:
    std::vector<float> dist;
    std::vector<int> parent;
    std::vector<bool> affected;
    std::vector<bool> affected_del;

    // OpenCL data structures
    bool opencl_available = false;
    OpenCLContext opencl_ctx;
    std::vector<std::pair<int, int>> edge_pairs;
    std::vector<float> edge_weights;

    SSSP(int V);
    void initialize(int source);
    void updateStep1(const Graph &graph, const std::vector<Edge> &inserts,
                     const std::vector<Edge> &deletes, bool use_openmp);
    void updateStep2(Graph &graph, bool use_openmp, int async_level, bool use_opencl = false);
    void updateStep2CPU(Graph &graph, bool use_openmp, int async_level); // Added declaration
    bool hasConverged(MPI_Comm comm);
    void markAffectedSubtree(int root, Graph &graph);

    // New method to prepare graph data for OpenCL
    void prepareGraphForOpenCL(const Graph &graph);
};

#endif // SSSP_H