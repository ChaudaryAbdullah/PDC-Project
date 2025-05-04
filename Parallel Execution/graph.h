#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <mpi.h>

struct Edge
{
    int u, v;
    float weight;
};

class Graph
{
public:
    int V, E;
    std::vector<Edge> edges;
    std::vector<std::vector<std::pair<int, float>>> adj;

    // Partitioning information
    std::vector<int> part;
    std::vector<int> local_vertices;
    std::vector<int> ghost_vertices;

    Graph();
    void loadFromFile(const std::string &filename);
    void partitionGraph(int num_parts);
    void distributeGraph(MPI_Comm comm);
    void addEdge(int u, int v, float weight);
    void applyUpdates(const std::vector<Edge> &updates);
    void gatherSSSPResults(MPI_Comm comm, std::vector<float> &global_dist);
};

#endif // GRAPH_H