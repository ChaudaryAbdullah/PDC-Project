#include "sssp.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <mpi.h>
#include <omp.h>
#include <iostream>

SSSP::SSSP(int V) : dist(V, std::numeric_limits<float>::infinity()),
                    parent(V, -1),
                    affected(V, false),
                    affected_del(V, false) {}

void SSSP::initialize(int source)
{
    if (source < 0 || source >= static_cast<int>(dist.size()))
    {
        std::cerr << "Error: Invalid source vertex " << source << std::endl;
        return;
    }

    std::fill(dist.begin(), dist.end(), std::numeric_limits<float>::infinity());
    std::fill(parent.begin(), parent.end(), -1);
    std::fill(affected.begin(), affected.end(), false);
    std::fill(affected_del.begin(), affected_del.end(), false);
    dist[source] = 0;

    // Mark source as affected to trigger initial computation
    affected[source] = true;
}

void SSSP::prepareGraphForOpenCL(const Graph &graph)
{
    // Convert adjacency list to edge pairs and weights for OpenCL
    edge_pairs.clear();
    edge_weights.clear();

    for (int u = 0; u < graph.V; u++)
    {
        for (const auto &neighbor : graph.adj[u])
        {
            int v = neighbor.first;
            float weight = neighbor.second;

            // Avoid duplicate edges (since this is an undirected graph)
            if (u < v)
            {
                edge_pairs.push_back({u, v});
                edge_weights.push_back(weight);
            }
        }
    }

    // Initialize OpenCL if not already done
    if (!opencl_available)
    {
        opencl_available = setupOpenCL(opencl_ctx, "relax_edges.cl");
        if (!opencl_available)
        {
            std::cerr << "Warning: OpenCL initialization failed, falling back to CPU implementation" << std::endl;
        }
    }
}

void SSSP::updateStep1(const Graph &graph, const std::vector<Edge> &inserts,
                       const std::vector<Edge> &deletes, bool use_openmp)
{
#pragma omp parallel for if (use_openmp)
    for (size_t i = 0; i < deletes.size(); i++)
    {
        const Edge &e = deletes[i];
        if (e.u < 0 || e.u >= static_cast<int>(dist.size()) ||
            e.v < 0 || e.v >= static_cast<int>(dist.size()))
        {
#pragma omp critical
            {
                std::cerr << "Warning: Invalid edge in deletions: " << e.u << " " << e.v << std::endl;
            }
            continue;
        }
        // Always mark both vertices as affected for deletion
        affected_del[e.u] = true;
        affected_del[e.v] = true;
        affected[e.u] = true;
        affected[e.v] = true;
        // Reset distances if the edge was part of the shortest path
        if (parent[e.v] == e.u)
        {
            dist[e.v] = std::numeric_limits<float>::infinity();
            parent[e.v] = -1;
        }
        if (parent[e.u] == e.v)
        {
            dist[e.u] = std::numeric_limits<float>::infinity();
            parent[e.u] = -1;
        }
    }

#pragma omp parallel for if (use_openmp)
    for (size_t i = 0; i < inserts.size(); i++)
    {
        const Edge &e = inserts[i];
        int u = e.u, v = e.v;
        float weight = e.weight;

        // Validate edge vertices
        if (u < 0 || u >= static_cast<int>(dist.size()) ||
            v < 0 || v >= static_cast<int>(dist.size()))
        {
#pragma omp critical
            {
                std::cerr << "Warning: Invalid edge in insertions: " << u << " " << v << std::endl;
            }
            continue;
        }

        if (dist[u] > dist[v])
            std::swap(u, v);

        if (dist[v] > dist[u] + weight)
        {
            dist[v] = dist[u] + weight;
            parent[v] = u;
            affected[v] = true;
        }
    }
}

void SSSP::updateStep2(Graph &graph, bool use_openmp, int async_level, bool use_opencl)
{
    if (use_opencl && opencl_available)
    {
        std::cout << "Running OpenCL SSSP on GPU..." << std::endl;
        prepareGraphForOpenCL(graph);
        runRelaxationKernel(opencl_ctx, dist, parent, edge_pairs, edge_weights);
    }
    else
    {
        std::cout << "Running CPU SSSP..." << std::endl;
        updateStep2CPU(graph, use_openmp, async_level);
    }
}

void SSSP::updateStep2CPU(Graph &graph, bool use_openmp, int async_level)
{
    bool changed;
    int iterations = 0;
    const int MAX_ITERATIONS = 100;

    // Preserve initial distances and reset only affected vertices
    std::vector<float> temp_dist = dist;
    std::vector<int> temp_parent = parent;

    do
    {
        changed = false;
        iterations++;

        // Phase 1: Handle deletions and reset affected subtrees
#pragma omp parallel for if (use_openmp) schedule(dynamic)
        for (int v = 0; v < graph.V; v++)
        {
            if (affected_del[v])
            {
                affected_del[v] = false;
                bool local_changed = false;
                for (const auto &neighbor : graph.adj[v])
                {
                    int c = neighbor.first;
                    if (c >= 0 && c < static_cast<int>(parent.size()) && parent[c] == v)
                    {
#pragma omp critical
                        {
                            dist[c] = std::numeric_limits<float>::infinity();
                            parent[c] = -1;
                            affected_del[c] = true;
                            affected[c] = true;
                            local_changed = true;
                        }
                    }
                }
                if (local_changed)
                {
#pragma omp atomic write
                    changed = true;
                }
            }
        }

        // Phase 2: Recompute paths for all vertices
        std::vector<bool> visited(graph.V, false);
        std::priority_queue<std::pair<float, int>, std::vector<std::pair<float, int>>, std::greater<>> pq;

        // Initialize with vertices that have finite distances
        for (int v = 0; v < graph.V; v++)
        {
            if (dist[v] != std::numeric_limits<float>::infinity())
            {
                pq.push({dist[v], v});
            }
        }

        while (!pq.empty())
        {
            auto [d, u] = pq.top();
            pq.pop();
            if (visited[u])
                continue;
            visited[u] = true;

            for (const auto &neighbor : graph.adj[u])
            {
                int v = neighbor.first;
                float weight = neighbor.second;
                if (v >= 0 && v < static_cast<int>(dist.size()))
                {
                    float new_dist = dist[u] + weight;
                    if (new_dist < dist[v])
                    {
                        dist[v] = new_dist;
                        parent[v] = u;
                        pq.push({new_dist, v});
                        affected[v] = true;
                        changed = true;
                    }
                }
            }
        }

        if (iterations % 10 == 0)
        {
            std::cout << "Iteration " << iterations << ", changed = " << changed << std::endl;
        }

        if (async_level > 1 && iterations % async_level == 0)
        {
            std::vector<float> global_dist(graph.V);
            for (int i = 0; i < graph.V; i++)
            {
                global_dist[i] = dist[i];
            }
            graph.gatherSSSPResults(MPI_COMM_WORLD, global_dist);

            for (int i = 0; i < graph.V; i++)
            {
                if (global_dist[i] < dist[i])
                {
                    dist[i] = global_dist[i];
                    affected[i] = true;
                    changed = true;
                }
            }
        }

    } while (changed && iterations < MAX_ITERATIONS);

    if (iterations >= MAX_ITERATIONS)
    {
        std::cerr << "Warning: updateStep2 reached maximum iterations without converging" << std::endl;
    }
    else
    {
        std::cout << "SSSP converged after " << iterations << " iterations." << std::endl;
    }
}

bool SSSP::hasConverged(MPI_Comm comm)
{
    int local_changed = 0;

    for (bool a : affected)
    {
        if (a)
        {
            local_changed = 1;
            break;
        }
    }

    int global_changed;
    MPI_Allreduce(&local_changed, &global_changed, 1, MPI_INT, MPI_LOR, comm);
    return global_changed == 0;
}

void SSSP::markAffectedSubtree(int root, Graph &graph)
{
    if (root < 0 || root >= static_cast<int>(dist.size()))
    {
        std::cerr << "Error: Invalid root vertex " << root << std::endl;
        return;
    }

    std::queue<int> q;
    q.push(root);
    affected_del[root] = true;
    affected[root] = true;

    while (!q.empty())
    {
        int v = q.front();
        q.pop();

        for (const auto &neighbor : graph.adj[v])
        {
            int c = neighbor.first;

            if (c < 0 || c >= static_cast<int>(parent.size()))
            {
                continue;
            }

            if (parent[c] == v)
            {
                dist[c] = std::numeric_limits<float>::infinity();
                parent[c] = -1;
                affected_del[c] = true;
                affected[c] = true;
                q.push(c);
            }
        }
    }
}