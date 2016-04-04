#include <tclap/CmdLine.h>

#include <andres/graphics/graphics-hdf5.hxx>

#include "lineage/graphics.hxx"

struct Parameters {
    std::string edgesFileName;
    std::string nodesFileName;
    std::string fragmentEdgeLabelsFileName;
    std::string outputFileName;
};

Parameters parseCommandLine(int argc, char** argv)
try
{
    Parameters parameters = {};

    TCLAP::CmdLine tclap("draw", ' ', "1.0");
    TCLAP::ValueArg<std::string> argNodesFileName("n", "nodes-file", "nodes information", true, parameters.nodesFileName, "nodes-file", tclap);
    TCLAP::ValueArg<std::string> argEdgesFileName("e", "edges-file", "edges information", true, parameters.edgesFileName, "edges-file", tclap);
    TCLAP::ValueArg<std::string> argFragmentEdgeLabelsFileName("s", "fragment-edge-labels-file", "solution", false, parameters.fragmentEdgeLabelsFileName, "fragment-edge-labels-file", tclap);
    TCLAP::ValueArg<std::string> argOutputFileName("o", "output-file", "output file", true, parameters.outputFileName, "output-hdf5-file", tclap);

    tclap.parse(argc, argv);

    parameters.edgesFileName = argEdgesFileName.getValue();
    parameters.nodesFileName = argNodesFileName.getValue();
    parameters.fragmentEdgeLabelsFileName = argFragmentEdgeLabelsFileName.getValue();
    parameters.outputFileName = argOutputFileName.getValue();

    return parameters;
}
catch(TCLAP::ArgException& e)
{
    throw std::runtime_error(e.error());
}

int main(int argc, char** argv)
try
{
    auto parameters = parseCommandLine(argc, argv);

    // load problem:
    auto problem = lineage::loadProblem(parameters.nodesFileName, parameters.edgesFileName);
    lineage::ProblemGraph problemGraph(problem);

    // draw:
    andres::graphics::Graphics<> graphics;
    if(parameters.fragmentEdgeLabelsFileName.empty()) {
        lineage::draw(problemGraph, graphics);
    }
    else {
        auto solution = lineage::loadSolution(parameters.fragmentEdgeLabelsFileName);
        lineage::SolutionGraph solutionGraph(problemGraph, solution);
        lineage::draw(solutionGraph, graphics);
    }

    // save graphics:
    hid_t file = andres::hdf5::createFile(parameters.outputFileName);
    andres::hdf5::save(file, graphics);
    andres::hdf5::closeFile(file);

    return 0;
}
catch (const std::runtime_error& error)
{
    std::cerr << "error: " << error.what() << std::endl;
    return 1;
}
