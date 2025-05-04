#ifndef UTILS_H
#define UTILS_H
#include "graph.h"
#include <vector>
#include <string>

std::vector<Edge> loadUpdates(const std::string &filename);
void saveResults(const std::string &filename, const std::vector<float> &dist);
void printStats(const std::vector<float> &dist);

#endif // UTILS_H