#pragma once
#ifndef LINEAGE_PROBLEM_HXX
#define LINEAGE_PROBLEM_HXX

#include <stdexcept>
#include <cmath>
#include <vector>
#include <fstream>


namespace lineage {

struct Node
{
    int t;
    int id;
    int cx, cy;
    double probability_birth_termination;
};

struct Edge
{
    int t0, v0;
    int t1, v1;
    double weight;
};

struct Problem
{
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<size_t> node_offsets;
};

inline
size_t loadNodes(const std::string& fileName, Problem& problem)
{
    problem.nodes.clear();
    problem.node_offsets.clear();
    problem.node_offsets.push_back(0);

    std::ifstream file(fileName);

    Node node;
    size_t counter = 0;
    while (file >> node.t >> node.id >> node.cx >> node.cy >> node.probability_birth_termination)
    {
        if (node.t == problem.node_offsets.size())
            problem.node_offsets.push_back(counter);

        problem.nodes.push_back(node);
        ++counter;
    }

    problem.node_offsets.push_back(counter);

    file.close();

    return counter;
}

inline
size_t loadEdges(const std::string& fileName, Problem& problem)
{
    problem.edges.clear();

    std::ifstream file(fileName);

    Edge e;
    size_t counter = 0;
    while (file >> e.t0 >> e.v0 >> e.t1 >> e.v1 >> e.weight)
    {
        problem.edges.push_back(e);
        ++counter;
    }

    file.close();

    return counter;
}

inline
Problem loadProblem(const std::string& nodesFileName, const std::string& edgesFileName)
{
    Problem problem;

    loadNodes(nodesFileName, problem);
    
    loadEdges(edgesFileName, problem);
    
    return problem;
}

template<class T = double>
struct NegativeLogProbabilityRatio
{
    typedef T value_type;

    NegativeLogProbabilityRatio(value_type epsilon = static_cast<value_type>(1) / static_cast<value_type>(255))
        :   epsilon_(epsilon),
            oneMinusEpsilon_(1.0 - epsilon)
        {
            if(epsilon <= .0 || epsilon * 2.0 >= 1.0) {
                throw std::out_of_range("epsilon out of range (0, 0.5).");
            }
        }
    value_type operator()(value_type x) const
        {
            assert(.0 <= x && x <= 1.0);
            if(x < epsilon_) {
                x = epsilon_;
            }
            else if(x > oneMinusEpsilon_) {
                x = oneMinusEpsilon_;
            }
            return std::log( (1.0 - x) / x );
        }

private:
    value_type epsilon_;
    value_type oneMinusEpsilon_;
};

} // namespace lineage

#endif
