#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include <deque>

#include <tclap/CmdLine.h>

#include "lineage/solution.hxx"
#include "lineage/solution-graph.hxx"

using namespace std;


struct Parameters {
    string edgesFileName;
    string nodesFileName;
    string solutionName;
};

Parameters parseCommandLine(int argc, char** argv)
try
{
    Parameters parameters;

    TCLAP::CmdLine tclap("track", ' ', "1.0");
    TCLAP::ValueArg<string> argNodesFileName("n", "nodes-file", "nodes information", true, parameters.nodesFileName, "nodes-file", tclap);
    TCLAP::ValueArg<string> argEdgesFileName("e", "edges-file", "edges information", true, parameters.edgesFileName, "edges-file", tclap);
    TCLAP::ValueArg<string> argSolutionName("s", "solution-name", "solution name", true, parameters.solutionName, "solution-name", tclap);
    
    tclap.parse(argc, argv);

    parameters.edgesFileName = argEdgesFileName.getValue();
    parameters.nodesFileName = argNodesFileName.getValue();
    parameters.solutionName = argSolutionName.getValue();

    return parameters;
}
catch(TCLAP::ArgException& e)
{
    throw runtime_error(e.error());
}

template<typename GRAPH, typename VLA, typename ELA>
inline void
edgeToVertexLabels(
    const GRAPH& graph,
    const ELA& edge_labels,
    VLA& vertex_labels
) {
    struct mask
    {
        mask(const ELA& edge_labels) : edge_labels_(edge_labels)
            {}
        bool vertex(std::size_t i) const
            { return true; }
        bool edge(std::size_t i) const
            { return !edge_labels_[i]; }

        const ELA& edge_labels_;
    };

    andres::graph::DepthFirstSearchData<> dfs_data(graph.numberOfVertices());
    for (std::size_t i = 0, label = 0; i < graph.numberOfVertices(); ++i)
        if (!dfs_data.visited(i))
        {
            andres::graph::depthFirstSearch(
                graph,
                mask(edge_labels),
                i,
                [&](std::size_t v, bool& proceed, bool& add_neighbors)
                {
                    vertex_labels[v] = label;
                    proceed = true;
                    add_neighbors = true;
                },
                dfs_data);

            ++label;
        }
}

int main(int argc, char** argv)
try
{
    auto parameters = parseCommandLine(argc, argv);

    // load problem:
    auto problem = lineage::loadProblem(parameters.nodesFileName, parameters.edgesFileName);

    lineage::ProblemGraph problemGraph(problem);

    auto solution = lineage::loadSolution(parameters.solutionName);
    
    lineage::SolutionGraph solutionGraph(problemGraph, solution);

    parameters.solutionName = parameters.solutionName.substr(0, parameters.solutionName.size() - 4);
    solutionGraph.save(parameters.solutionName);
    solutionGraph.saveSVG(parameters.solutionName + "-lineage-tree.svg");

    return 0;
}
catch(const runtime_error& error) {
    cerr << "error: " << error.what() << endl;
    return 1;
}
