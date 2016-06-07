#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include <deque>

#include <tclap/CmdLine.h>

#include "lineage/validation.hxx"

using namespace std;

struct Parameters {
    string edgesFileName;
    string nodesFileName;
    string fragmentEdgeLabelsFileName;
};

Parameters parseCommandLine(int argc, char** argv)
try
{
    Parameters parameters;

    TCLAP::CmdLine tclap("track", ' ', "1.0");
    TCLAP::ValueArg<string> argNodesFileName("n", "nodes-file", "nodes information", true, parameters.nodesFileName, "nodes-file", tclap);
    TCLAP::ValueArg<string> argEdgesFileName("e", "edges-file", "edges information", true, parameters.edgesFileName, "edges-file", tclap);
    TCLAP::ValueArg<string> argFragmentEdgeLabelsFileName("s", "fragment-edge-labels-file", "solution", false, parameters.fragmentEdgeLabelsFileName, "fragment-edge-labels-file", tclap);

    tclap.parse(argc, argv);

    parameters.edgesFileName = argEdgesFileName.getValue();
    parameters.nodesFileName = argNodesFileName.getValue();
    parameters.fragmentEdgeLabelsFileName = argFragmentEdgeLabelsFileName.getValue();

    return parameters;
}
catch (TCLAP::ArgException& e)
{
    throw runtime_error(e.error());
}

int main(int argc, char** argv)
try
{
    auto parameters = parseCommandLine(argc, argv);
    auto problem = lineage::loadProblem(parameters.nodesFileName, parameters.edgesFileName);
    auto solution = lineage::loadSolution(parameters.fragmentEdgeLabelsFileName);

    lineage::ProblemGraph problemGraph(problem);
    lineage::validate(problemGraph, solution);

    return 0;
}
catch (const runtime_error& error)
{
    cerr << "error: " << error.what() << endl;
    return 1;
}
