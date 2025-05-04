#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>
#include <limits>
#include "graph.h"
#include "sssp.h"
#include "utils.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4)
    {
        if (rank == 0)
        {
            std::cerr << "Usage: " << argv[0]
                      << " <graph_file> <updates_file> <source_vertex> [output_file] [--openmp] [--async=<level>] [--opencl]" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    std::string graph_file = argv[1];
    std::string updates_file = argv[2];

    // Check if the source vertex is a valid number
    int source;
    try
    {
        source = std::stoi(argv[3]);
    }
    catch (const std::exception &e)
    {
        if (rank == 0)
        {
            std::cerr << "Error: Source vertex must be a valid integer, got '" << argv[3] << "'" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Default values
    std::string output_file = "";
    bool use_openmp = false;
    bool use_opencl = false;
    int async_level = 1;

    // Process optional arguments
    for (int i = 4; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "--openmp")
        {
            use_openmp = true;
        }
        else if (arg == "--opencl")
        {
            use_opencl = true;
        }
        else if (arg.compare(0, 8, "--async=") == 0)
        {
            try
            {
                async_level = std::stoi(arg.substr(8));
                if (async_level <= 0)
                {
                    if (rank == 0)
                    {
                        std::cerr << "Warning: Invalid async level " << async_level << ", setting to 1" << std::endl;
                    }
                    async_level = 1;
                }
            }
            catch (const std::exception &e)
            {
                if (rank == 0)
                {
                    std::cerr << "Warning: Invalid async value, using default level 1" << std::endl;
                }
                async_level = 1;
            }
        }
        else if (arg.compare(0, 2, "--") != 0)
        {
            output_file = arg;
        }
        else
        {
            if (rank == 0)
            {
                std::cerr << "Warning: Unknown option '" << arg << "'" << std::endl;
            }
        }
    }

    if (rank == 0)
    {
        std::cout << "Configuration:" << std::endl;
        std::cout << "  Graph file: " << graph_file << std::endl;
        std::cout << "  Updates file: " << updates_file << std::endl;
        std::cout << "  Source vertex: " << source << std::endl;
        std::cout << "  Output file: " << (output_file.empty() ? "none" : output_file) << std::endl;
        std::cout << "  OpenMP: " << (use_openmp ? "enabled" : "disabled") << std::endl;
        std::cout << "  OpenCL: " << (use_opencl ? "enabled" : "disabled") << std::endl;
        std::cout << "  Async level: " << async_level << std::endl;
    }

    // Load graph and partition it
    Graph graph;
    if (rank == 0)
    {
        std::cout << "Loading graph from " << graph_file << std::endl;
        graph.loadFromFile(graph_file);
        std::cout << "Graph loaded: " << graph.V << " vertices, " << graph.E << " edges" << std::endl;
        graph.partitionGraph(size);
    }

    // First broadcast graph size to all processes
    int graph_info[2];
    if (rank == 0)
    {
        graph_info[0] = graph.V;
        graph_info[1] = graph.E;
    }
    MPI_Bcast(graph_info, 2, MPI_INT, 0, MPI_COMM_WORLD);

    // Now non-root processes know the graph size
    if (rank != 0)
    {
        graph.V = graph_info[0];
        graph.E = graph_info[1];
        graph.adj.resize(graph.V);
        graph.part.resize(graph.V);
    }

    // Broadcast partition information
    MPI_Bcast(graph.part.data(), graph.V, MPI_INT, 0, MPI_COMM_WORLD);

    // Distribute graph data
    std::vector<Edge> all_edges;
    if (rank == 0)
    {
        all_edges = graph.edges;
    }

    int num_edges = (rank == 0) ? graph.E : 0;
    MPI_Bcast(&num_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        all_edges.resize(num_edges);
    }

    // Create custom MPI datatype for Edge
    MPI_Datatype MPI_EDGE;
    int blocklengths[3] = {1, 1, 1};
    MPI_Aint offsets[3];
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_FLOAT};

    Edge temp;
    MPI_Aint base_address;
    MPI_Get_address(&temp, &base_address);
    MPI_Get_address(&temp.u, &offsets[0]);
    MPI_Get_address(&temp.v, &offsets[1]);
    MPI_Get_address(&temp.weight, &offsets[2]);
    offsets[0] = MPI_Aint_diff(offsets[0], base_address);
    offsets[1] = MPI_Aint_diff(offsets[1], base_address);
    offsets[2] = MPI_Aint_diff(offsets[2], base_address);

    MPI_Type_create_struct(3, blocklengths, offsets, types, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);

    // Broadcast edges
    MPI_Bcast(all_edges.data(), num_edges, MPI_EDGE, 0, MPI_COMM_WORLD);

    // Reconstruct graph if not root
    if (rank != 0)
    {
        for (const auto &edge : all_edges)
        {
            graph.adj[edge.u].emplace_back(edge.v, edge.weight);
            graph.adj[edge.v].emplace_back(edge.u, edge.weight);
        }
        graph.edges = all_edges;
    }

    graph.distributeGraph(MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "Graph distributed. Process 0 has " << graph.local_vertices.size()
                  << " local vertices and " << graph.ghost_vertices.size() << " ghost vertices" << std::endl;
    }

    // Initialize SSSP
    SSSP sssp(graph.V);
    sssp.initialize(source);

    if (rank == 0)
    {
        std::cout << "Running initial SSSP calculation from source " << source << std::endl;
    }

    sssp.updateStep2(graph, use_openmp, async_level, use_opencl);

    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<float> initial_dist(graph.V);
    for (int i = 0; i < graph.V; i++)
    {
        initial_dist[i] = sssp.dist[i];
    }

    graph.gatherSSSPResults(MPI_COMM_WORLD, initial_dist);

    if (rank == 0)
    {
        std::cout << "Initial SSSP completed. Statistics:" << std::endl;
        printStats(initial_dist);
    }

    // Load updates
    std::vector<Edge> all_updates;
    if (rank == 0)
    {
        std::cout << "Loading updates from " << updates_file << std::endl;
        all_updates = loadUpdates(updates_file);
        std::cout << "Loaded " << all_updates.size() << " updates" << std::endl;
    }

    int num_updates = (rank == 0) ? all_updates.size() : 0;
    MPI_Bcast(&num_updates, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        all_updates.resize(num_updates);
    }

    MPI_Bcast(all_updates.data(), num_updates, MPI_EDGE, 0, MPI_COMM_WORLD);

    std::vector<Edge> inserts, deletes;
    for (const auto &e : all_updates)
    {
        if (e.weight >= 0)
        {
            inserts.push_back(e);
        }
        else
        {
            Edge delete_edge = {e.u, e.v, -1.0f};
            for (const auto &adj_edge : graph.adj[e.u])
            {
                if (adj_edge.first == e.v)
                {
                    delete_edge.weight = adj_edge.second;
                    break;
                }
            }
            deletes.push_back(delete_edge);
        }
    }

    if (rank == 0)
    {
        std::cout << "Processing " << inserts.size() << " insertions and "
                  << deletes.size() << " deletions" << std::endl;
    }

    double start_time = MPI_Wtime();

    graph.applyUpdates(all_updates); // Apply both insertions and deletions

    // Redistribute updated graph
    num_edges = graph.E;
    std::vector<Edge> updated_edges;
    if (rank == 0)
    {
        updated_edges.clear();
        for (int u = 0; u < graph.V; u++)
        {
            for (const auto &neighbor : graph.adj[u])
            {
                int v = neighbor.first;
                float weight = neighbor.second;
                if (u < v) // Avoid duplicates in undirected graph
                {
                    updated_edges.push_back({u, v, weight});
                }
            }
        }
        num_edges = updated_edges.size();
    }
    MPI_Bcast(&num_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        updated_edges.resize(num_edges);
    }
    MPI_Bcast(updated_edges.data(), num_edges, MPI_EDGE, 0, MPI_COMM_WORLD);
    if (rank != 0)
    {
        graph.adj.clear();
        graph.adj.resize(graph.V);
        for (const auto &edge : updated_edges)
        {
            graph.adj[edge.u].emplace_back(edge.v, edge.weight);
            graph.adj[edge.v].emplace_back(edge.u, edge.weight);
        }
        graph.edges = updated_edges;
        graph.E = num_edges;
    }
    graph.distributeGraph(MPI_COMM_WORLD);

    sssp.updateStep1(graph, inserts, deletes, use_openmp);

    MPI_Barrier(MPI_COMM_WORLD);

    sssp.updateStep2(graph, use_openmp, async_level, use_opencl);

    std::vector<float> global_dist(graph.V, std::numeric_limits<float>::infinity());
    for (int i = 0; i < graph.V; i++)
    {
        global_dist[i] = sssp.dist[i];
    }
    graph.gatherSSSPResults(MPI_COMM_WORLD, global_dist);

    double end_time = MPI_Wtime();

    if (rank == 0)
    {
        std::cout << "SSSP update completed in " << (end_time - start_time) << " seconds\n";
        printStats(global_dist);

        if (!output_file.empty())
        {
            saveResults(output_file, global_dist);
            std::cout << "Results saved to " << output_file << "\n";
        }
    }

    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    return 0;
}