#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>

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
        file << i << " " << std::fixed << std::setprecision(2) << dist[i] << "\n";
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