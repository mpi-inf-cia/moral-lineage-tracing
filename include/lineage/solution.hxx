#pragma once
#ifndef LINEAGE_SOLUTION_HXX
#define LINEAGE_SOLUTION_HXX

#include <fstream>

#include <andres/graph/dfs.hxx>

namespace lineage {

struct Solution {
    typedef std::vector<unsigned char> EdgeLabels;

    EdgeLabels edge_labels;
};

inline
void saveSolution(std::string const& fileName, Solution const& solution)
{
    std::ofstream file(fileName);

    for (size_t j = 0; j < solution.edge_labels.size(); ++j)
        file << static_cast<size_t>(solution.edge_labels[j]) << std::endl;
    
    file.close();
}

inline
Solution loadSolution(std::string const& fileName)
{
    Solution solution;

    std::ifstream file(fileName);
    
    size_t edgeLabel = 0;
    while (file >> edgeLabel)
        solution.edge_labels.push_back(edgeLabel);
    
    file.close();

    return solution;
}

} // namespace lineage

#endif
