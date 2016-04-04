#include <stdexcept>

#include <andres/ilp/gurobi-callback.hxx>
// #include <andres/ilp/gurobi.hxx>

#include <tclap/CmdLine.h>

#include "lineage/solver-ilp.hxx"
#include "lineage/solution-graph.hxx"

struct Parameters {
    std::string edgesFileName;
    std::string nodesFileName;
    std::string solutionName;
    double biasSpatial { .5 };
    double biasTemporal { .5 };
    double terminationCost { .0 };
    double birthCost { .0 };
    bool bifurcationConstraint { false };
};

Parameters parseCommandLine(int argc, char** argv)
try
{
    Parameters parameters;

    TCLAP::CmdLine tclap("track", ' ', "1.0");
    TCLAP::ValueArg<std::string> argNodesFileName("n", "nodes-file", "nodes information", true, parameters.nodesFileName, "nodes-file", tclap);
    TCLAP::ValueArg<std::string> argEdgesFileName("e", "edges-file", "edges information", true, parameters.edgesFileName, "edges-file", tclap);
    TCLAP::ValueArg<std::string> argSolutionName("s", "solution-name", "solution name", true, parameters.solutionName, "solution-name", tclap);
    TCLAP::ValueArg<double> argBiasSpatial("b", "cut-prior-spatial", "cut prior spatial", false, parameters.biasSpatial, "cut prior spatial", tclap);
    TCLAP::ValueArg<double> argBiasTemporal("t", "cut-prior-temporal", "cut prior temporal", false, parameters.biasTemporal, "cut prior temporal", tclap);
    TCLAP::ValueArg<double> argTerminationCost("T", "termination-cost", "early termination cost", false, parameters.terminationCost, "early termination cost", tclap);
    TCLAP::ValueArg<double> argBirthCost("B", "birth-cost", "birth cost", false, parameters.birthCost, "birth cost", tclap);
    TCLAP::SwitchArg argBifurcationConstraint("F", "bifurcation-constraint", "Enforce bifurcation contraint. (Default: disabled).", tclap);
    
    tclap.parse(argc, argv);

    parameters.edgesFileName = argEdgesFileName.getValue();
    parameters.nodesFileName = argNodesFileName.getValue();
    parameters.solutionName = argSolutionName.getValue();
    parameters.biasSpatial = argBiasSpatial.getValue();
    parameters.biasTemporal = argBiasTemporal.getValue();
    parameters.terminationCost = argTerminationCost.getValue();
    parameters.birthCost = argBirthCost.getValue();
    parameters.bifurcationConstraint = argBifurcationConstraint.getValue();

    if (parameters.biasSpatial < std::numeric_limits<double>::epsilon() || parameters.biasSpatial > 1.0 - std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Spatial bias must be in the range (0, 1)");

    if (parameters.biasTemporal < std::numeric_limits<double>::epsilon() || parameters.biasTemporal > 1.0 - std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Temporal bias must be in the range (0, 1)");

    std::cout << "parameters:" << std::endl
        << "- cut prior (spatial) : " << parameters.biasSpatial << std::endl
        << "- cut prior (temporal): " << parameters.biasTemporal << std::endl
        << "- cost of termination : " << parameters.terminationCost << std::endl
        << "- cost of birth: " << parameters.birthCost << std::endl
        << "- bifurcation constraint: " << (parameters.bifurcationConstraint ? "yes" : "no") << std::endl
        << std::endl;

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

    // map edge probabilities to edge cut costs:
    lineage::NegativeLogProbabilityRatio<> func;
    for (auto& e : problem.edges)
        if (e.t0 != e.t1)
            e.weight = func(e.weight) + func(parameters.biasTemporal);
        else
            e.weight = func(e.weight) + func(parameters.biasSpatial);

    // solve problem:
    auto solution = lineage::solver_ilp<andres::ilp::Gurobi>(
        problem,
        parameters.terminationCost,
        parameters.birthCost,
        parameters.bifurcationConstraint,
        parameters.solutionName
    );
    
    // save solution:
    lineage::ProblemGraph problemGraph(problem);
    lineage::SolutionGraph solutionGraph(problemGraph, solution);
    solutionGraph.save(parameters.solutionName);
    solutionGraph.saveSVG(parameters.solutionName + "-lineage-tree.svg");

    return 0;
}
catch (const std::runtime_error& error)
{
    std::cerr << "error: " << error.what() << std::endl;
    return 1;
}
