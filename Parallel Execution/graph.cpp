#include "graph.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <metis.h>
#include <algorithm>

Graph::Graph() : V(0), E(0) {}

void Graph::loadFromFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    if (!(file >> V >> E))
    {
        std::cerr << "Error reading graph header" << std::endl;
        return;
    }

    if (V <= 0 || E <= 0)
    {
        std::cerr << "Invalid graph size: V=" << V << ", E=" << E << std::endl;
        V = 0;
        E = 0;
        return;
    }

    adj.clear();
    adj.resize(V);
    edges.clear();
    edges.reserve(E);

    int u, v;
    float weight;
    int edge_count = 0;
    std::string line;

    std::getline(file, line);

    while (edge_count < E && std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);

        if (!(iss >> u >> v >> weight))
        {
            std::cerr << "Error parsing edge line: " << line << std::endl;
            continue;
        }

        if (u < 0 || u >= V || v < 0 || v >= V)
        {
            std::cerr << "Invalid vertex indices in edge: " << u << " " << v << std::endl;
            continue;
        }

        if (u == v)
        {
            std::cerr << "Warning: Self-loop found at vertex " << u << ", ignoring" << std::endl;
            continue;
        }

        if (weight < 0)
        {
            std::cerr << "Warning: Negative weight found in edge " << u << "-" << v
                      << ", Dijkstra's algorithm may not work correctly" << std::endl;
        }

        edges.push_back({u, v, weight});
        adj[u].emplace_back(v, weight);
        adj[v].emplace_back(u, weight);
        edge_count++;
    }

    if (edge_count < E)
    {
        std::cerr << "Warning: Expected " << E << " edges but found only " << edge_count << std::endl;
        E = edge_count;
    }

    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);

        if (iss >> u >> v >> weight)
        {
            if (u < 0 || u >= V || v < 0 || v >= V)
            {
                std::cerr << "Invalid vertex indices in additional edge: " << u << " " << v << std::endl;
                continue;
            }

            if (u == v)
            {
                std::cerr << "Warning: Self-loop found at vertex " << u << ", ignoring" << std::endl;
                continue;
            }

            edges.push_back({u, v, weight});
            adj[u].emplace_back(v, weight);
            adj[v].emplace_back(u, weight);
            E++;
        }
    }

    std::cout << "Successfully loaded graph with " << V << " vertices and " << E << " edges" << std::endl;
}

void Graph::partitionGraph(int num_parts)
{
    if (V == 0)
        return;
    if (num_parts > V)
    {
        std::cerr << "Warning: More partitions than vertices, setting num_parts = V" << std::endl;
        num_parts = V;
    }

    idx_t nvtxs = V;
    idx_t ncon = 1;
    idx_t *xadj = new idx_t[V + 1];
    idx_t *adjncy = new idx_t[2 * E];
    idx_t *vwgt = nullptr;
    idx_t *adjwgt = nullptr;
    idx_t objval;
    idx_t nparts = num_parts;
    part.resize(V);

    xadj[0] = 0;
    for (int i = 0; i < V; i++)
    {
        xadj[i + 1] = xadj[i] + adj[i].size();
        for (size_t j = 0; j < adj[i].size(); j++)
        {
            adjncy[xadj[i] + j] = adj[i][j].first;
        }
    }

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    
    // Remove the contiguity constraint since the graph might not be contiguous
    // options[METIS_OPTION_CONTIG] = 1;  // This was causing the error
    
    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, vwgt, nullptr,
                                  adjwgt, &nparts, nullptr, nullptr, options,
                                  &objval, part.data());

    if (ret != METIS_OK)
    {
        std::cerr << "METIS partitioning failed with code " << ret << std::endl;
        
        // If partitioning still fails, provide more diagnostics
        if (ret == METIS_ERROR_INPUT)
            std::cerr << "Error in the graph's input format" << std::endl;
        else if (ret == METIS_ERROR_MEMORY)
            std::cerr << "METIS could not allocate required memory" << std::endl;
        else if (ret == METIS_ERROR)
            std::cerr << "General METIS error" << std::endl;
            
        std::cerr << "Using simple vertex partitioning instead" << std::endl;
        for (int v = 0; v < V; v++)
        {
            part[v] = v % num_parts;
        }
    }
    else
    {
        std::cout << "Successfully partitioned graph into " << num_parts << " parts" << std::endl;
        
        // Optionally, print partition statistics
        std::vector<int> part_sizes(num_parts, 0);
        for (int v = 0; v < V; v++)
        {
            if (part[v] >= 0 && part[v] < num_parts)
                part_sizes[part[v]]++;
        }
        
        std::cout << "Partition sizes: ";
        for (int i = 0; i < num_parts; i++)
        {
            std::cout << part_sizes[i];
            if (i < num_parts - 1)
                std::cout << ", ";
        }
        std::cout << std::endl;
    }

    delete[] xadj;
    delete[] adjncy;
}
void Graph::distributeGraph(MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    local_vertices.clear();
    ghost_vertices.clear();

    for (int v = 0; v < V; v++)
    {
        if (part[v] == rank)
        {
            local_vertices.push_back(v);
        }
    }

    for (int v : local_vertices)
    {
        for (const auto &neighbor : adj[v])
        {
            int u = neighbor.first;
            if (part[u] != rank &&
                std::find(ghost_vertices.begin(), ghost_vertices.end(), u) == ghost_vertices.end())
            {
                ghost_vertices.push_back(u);
            }
        }
    }
}

void Graph::addEdge(int u, int v, float weight)
{
    if (u < 0 || u >= V || v < 0 || v >= V)
    {
        std::cerr << "Invalid vertex indices in edge: " << u << " " << v << std::endl;
        return;
    }

    edges.push_back({u, v, weight});
    adj[u].emplace_back(v, weight);
    adj[v].emplace_back(u, weight);
    E++;
}

void Graph::applyUpdates(const std::vector<Edge> &updates)
{
    for (const auto &edge : updates)
    {
        if (edge.weight < 0) // Handle deletion
        {
            // Remove edge from adj[edge.u]
            adj[edge.u].erase(
                std::remove_if(adj[edge.u].begin(), adj[edge.u].end(),
                               [edge](const std::pair<int, float> &neighbor)
                               {
                                   return neighbor.first == edge.v;
                               }),
                adj[edge.u].end());
            // Remove edge from adj[edge.v]
            adj[edge.v].erase(
                std::remove_if(adj[edge.v].begin(), adj[edge.v].end(),
                               [edge](const std::pair<int, float> &neighbor)
                               {
                                   return neighbor.first == edge.u;
                               }),
                adj[edge.v].end());
            // Update edge count
            E--;
        }
        else // Handle insertion or update
        {
            bool found = false;
            for (auto &neighbor : adj[edge.u])
            {
                if (neighbor.first == edge.v)
                {
                    neighbor.second = edge.weight;
                    found = true;
                    break;
                }
            }
            for (auto &neighbor : adj[edge.v])
            {
                if (neighbor.first == edge.u)
                {
                    neighbor.second = edge.weight;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                addEdge(edge.u, edge.v, edge.weight);
            }
        }
    }
}

void Graph::gatherSSSPResults(MPI_Comm comm, std::vector<float> &global_dist)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<float> local_dists = global_dist;
    std::vector<float> all_dists(V);

    MPI_Allreduce(local_dists.data(), all_dists.data(), V, MPI_FLOAT, MPI_MIN, comm);

    global_dist = all_dists;

    if (rank == 0)
    {
        std::cout << "Gathered SSSP results from all processes" << std::endl;
    }
}