#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <queue>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <tuple>

// Edge structure
struct Edge
{
    int u, v;
    float weight;
};

// Graph class
class Graph
{
public:
    int V, E;
    std::vector<Edge> edges;
    std::vector<std::vector<std::pair<int, float>>> adj;

    Graph() : V(0), E(0) {}

    void loadFromFile(const std::string &filename)
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

        std::getline(file, line); // Skip the first line after reading V and E

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
            adj[v].emplace_back(u, weight); // For undirected graph
            edge_count++;
        }

        if (edge_count < E)
        {
            std::cerr << "Warning: Expected " << E << " edges but found only " << edge_count << std::endl;
            E = edge_count;
        }

        std::cout << "Successfully loaded graph with " << V << " vertices and " << E << " edges" << std::endl;
    }

    void addEdge(int u, int v, float weight)
    {
        if (u < 0 || u >= V || v < 0 || v >= V)
        {
            std::cerr << "Invalid vertex indices in edge: " << u << " " << v << std::endl;
            return;
        }

        // Check if edge already exists in adj[u]
        for (auto &neighbor : adj[u])
        {
            if (neighbor.first == v)
            {
                neighbor.second = weight; // Update weight
                // Also update the reverse direction
                for (auto &rev_neighbor : adj[v])
                {
                    if (rev_neighbor.first == u)
                    {
                        rev_neighbor.second = weight;
                        return;
                    }
                }
                return;
            }
        }

        // Check if edge exists in adj[v] (reverse direction)
        for (auto &neighbor : adj[v])
        {
            if (neighbor.first == u)
            {
                neighbor.second = weight;
                // Add the forward direction since it wasn't found
                adj[u].emplace_back(v, weight);
                edges.push_back({u, v, weight});
                E++;
                return;
            }
        }

        // If edge doesn't exist in either direction, add it
        edges.push_back({u, v, weight});
        adj[u].emplace_back(v, weight);
        adj[v].emplace_back(u, weight); // For undirected graph
        E++;
    }

    void applyUpdates(const std::vector<Edge> &updates)
    {
        for (const auto &edge : updates)
        {
            if (edge.weight < 0)
            {
                // Handle edge deletion
                removeEdge(edge.u, edge.v);
                continue;
            }

            // Treat as insertion or update
            addEdge(edge.u, edge.v, edge.weight);
        }
    }

    void removeEdge(int u, int v)
    {
        if (u < 0 || u >= V || v < 0 || v >= V)
        {
            std::cerr << "Invalid vertex indices in edge deletion: " << u << " " << v << std::endl;
            return;
        }

        // Remove from edge list
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                   [u, v](const Edge &e)
                                   {
                                       return (e.u == u && e.v == v) || (e.u == v && e.v == u);
                                   }),
                    edges.end());

        // Remove from u's adjacency list
        adj[u].erase(std::remove_if(adj[u].begin(), adj[u].end(),
                                    [v](const std::pair<int, float> &p)
                                    {
                                        return p.first == v;
                                    }),
                     adj[u].end());

        // Remove from v's adjacency list
        adj[v].erase(std::remove_if(adj[v].begin(), adj[v].end(),
                                    [u](const std::pair<int, float> &p)
                                    {
                                        return p.first == u;
                                    }),
                     adj[v].end());

        E--;
    }
};

// SSSP (Single Source Shortest Path) class
class SSSP
{
public:
    std::vector<float> dist;
    std::vector<int> parent;

    SSSP(int V) : dist(V, std::numeric_limits<float>::infinity()),
                  parent(V, -1) {}

    void initialize(int source)
    {
        if (source < 0 || source >= static_cast<int>(dist.size()))
        {
            std::cerr << "Error: Invalid source vertex " << source << std::endl;
            return;
        }

        std::fill(dist.begin(), dist.end(), std::numeric_limits<float>::infinity());
        std::fill(parent.begin(), parent.end(), -1);
        dist[source] = 0;
    }

    void dijkstra(const Graph &graph, int source)
    {
        initialize(source);

        using P = std::pair<float, int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        pq.emplace(0, source);

        while (!pq.empty())
        {
            float d = pq.top().first;
            int u = pq.top().second;
            pq.pop();

            if (d > dist[u])
                continue;

            for (const auto &neighbor : graph.adj[u])
            {
                int v = neighbor.first;
                float weight = neighbor.second;

                if (dist[v] > dist[u] + weight)
                {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                    pq.emplace(dist[v], v);
                }
            }
        }

        std::cout << "Dijkstra's algorithm completed." << std::endl;
    }
};

// Utility functions
std::vector<Edge> loadUpdates(const std::string &filename)
{
    std::vector<Edge> updates;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening updates file: " << filename << std::endl;
        return updates;
    }

    std::string line;
    while (std::getline(file, line))
    {
        // Skip empty lines
        if (line.empty())
            continue;

        // Skip comment lines or lines that don't start with a digit
        if (line[0] == '#' || !isdigit(line[0]))
            continue;

        std::istringstream iss(line);
        Edge e;

        // Try to parse the line as "u v weight"
        if (iss >> e.u >> e.v)
        {
            std::string weight_str;
            if (iss >> weight_str)
            {
                // Check if it's a removal (marked with '-')
                if (weight_str == "-")
                {
                    e.weight = -1.0f; // Using -1 to indicate removal
                }
                else
                {
                    try
                    {
                        // Try to parse as float
                        e.weight = std::stof(weight_str);
                    }
                    catch (const std::exception &ex)
                    {
                        std::cerr << "Error parsing weight in line: " << line << std::endl;
                        continue; // Skip this malformed line
                    }
                }

                updates.push_back(e);
                std::cout << "Loaded update: " << e.u << " " << e.v << " " << e.weight << std::endl;
            }
        }
        else
        {
            std::cerr << "Malformed update line: " << line << std::endl;
        }
    }

    std::cout << "Total updates loaded: " << updates.size() << std::endl;
    return updates;
}

void saveResults(const std::string &filename, const std::vector<float> &dist)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening output file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < dist.size(); i++)
    {
        file << i << " " << std::fixed << std::setprecision(2);

        if (dist[i] == std::numeric_limits<float>::infinity())
            file << "inf\n";
        else
            file << dist[i] << "\n";
    }
}

void printStats(const std::vector<float> &dist)
{
    int reachable = 0;
    float max_dist = 0;
    float sum_dist = 0;

    for (float d : dist)
    {
        if (d < std::numeric_limits<float>::infinity())
        {
            reachable++;
            if (d > max_dist)
                max_dist = d;
            sum_dist += d;
        }
    }

    std::cout << "SSSP Statistics:\n";
    std::cout << "  Reachable vertices: " << reachable << "/" << dist.size() << "\n";
    std::cout << "  Maximum distance: " << max_dist << "\n";
    std::cout << "  Average distance: " << (reachable > 0 ? sum_dist / reachable : 0) << "\n";
}

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <graph_file> <updates_file> <source_vertex> [output_file]" << std::endl;
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
        std::cerr << "Error: Source vertex must be a valid integer, got '" << argv[3] << "'" << std::endl;
        return 1;
    }

    // Default values
    std::string output_file = "";

    // Process optional output file
    if (argc > 4)
    {
        output_file = argv[4];
    }

    std::cout << "Configuration:" << std::endl;
    std::cout << "  Graph file: " << graph_file << std::endl;
    std::cout << "  Updates file: " << updates_file << std::endl;
    std::cout << "  Source vertex: " << source << std::endl;
    std::cout << "  Output file: " << (output_file.empty() ? "none" : output_file) << std::endl;

    // Load graph
    Graph graph;
    std::cout << "Loading graph from " << graph_file << std::endl;
    graph.loadFromFile(graph_file);

    // Validate source vertex
    if (source < 0 || source >= graph.V)
    {
        std::cerr << "Error: Source vertex " << source << " is out of range (0 to " << (graph.V - 1) << ")" << std::endl;
        return 1;
    }

    std::cout << "Graph loaded: " << graph.V << " vertices, " << graph.E << " edges" << std::endl;

    // Initialize SSSP
    SSSP sssp(graph.V);

    std::cout << "Running initial SSSP calculation from source " << source << std::endl;
    sssp.dijkstra(graph, source);

    std::vector<float> initial_dist = sssp.dist;
    std::cout << "Initial SSSP completed. Statistics:" << std::endl;
    printStats(initial_dist);

    // Load updates
    std::cout << "Loading updates from " << updates_file << std::endl;
    std::vector<Edge> all_updates = loadUpdates(updates_file);
    std::cout << "Loaded " << all_updates.size() << " updates" << std::endl;

    std::cout << "Processing " << all_updates.size() << " updates" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    // Apply all updates to the graph structure
    graph.applyUpdates(all_updates);

    // Recompute SSSP from scratch
    sssp.dijkstra(graph, source);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

    std::cout << "SSSP update completed in " << duration << " seconds\n";
    printStats(sssp.dist);

    if (!output_file.empty())
    {
        saveResults(output_file, sssp.dist);
        std::cout << "Results saved to " << output_file << "\n";
    }

    return 0;
}