#pragma once
#ifndef LINEAGE_GRAPHICS_HXX
#define LINEAGE_GRAPHICS_HXX

#include <random>
#include <array>

#include <andres/graphics/graphics.hxx>

#include "solution-graph.hxx"

namespace lineage {

template<class T, class S>
void draw(ProblemGraph const& problemGraph, andres::graphics::Graphics<T, S>& graphics)
{
    for (size_t j = 0; j < problemGraph.problem().nodes.size(); ++j)
    {
        auto const& node = problemGraph.problem().nodes[j];

        graphics.definePoint(node.cx, node.cy, node.t);
    }

    for (size_t j = 0; j < problemGraph.problem().edges.size(); ++j)
    {
        auto const& edge = problemGraph.problem().edges[j];

        auto lineProperty = graphics.defineLineProperty(
            true,
            edge.weight * 255,
            (1.0 - edge.weight) * 255,
            0,
            (1.0 - edge.weight) * 255
        );

        graphics.defineLine(
            problemGraph.nodeInFrame(edge.t0, edge.v0),
            problemGraph.nodeInFrame(edge.t1, edge.v1),
            lineProperty
        );
    }

    graphics.center();
    graphics.normalize(0, 1);
    graphics.normalize(2);
}

template<class T, class S>
void draw(SolutionGraph const& solutionGraph, andres::graphics::Graphics<T, S>& graphics)
{
    Problem const& problem = solutionGraph.problem();
    ProblemGraph const& problemGraph = solutionGraph.problemGraph();
    Solution const& solution = solutionGraph.solution();

    // make as many random colors as there are cells
    typedef std::array<unsigned char, 3> Color;
    std::vector<Color> colors(solutionGraph.numberOfCells());

    int seed = 4711;
    std::default_random_engine randomEngine(seed);
    std::uniform_int_distribution<unsigned char> dist(0, 255);

    for (size_t j = 0; j < solutionGraph.numberOfCells(); ++j)
        for (size_t k = 0; k < 3; ++k)
            colors[j][k] = dist(randomEngine);

    for (size_t j = 0; j < problemGraph.problem().nodes.size(); ++j)
    {
        auto const& node = problemGraph.problem().nodes[j];
        
        graphics.definePoint(node.cx, node.cy, node.t);
    }

    for (size_t j = 0; j < problemGraph.problem().edges.size(); ++j)
    {
        auto const& edge = problemGraph.problem().edges[j];
        
        auto v0 = problemGraph.nodeInFrame(edge.t0, edge.v0);
        auto v1 = problemGraph.nodeInFrame(edge.t1, edge.v1);

        auto edge_label = solution.edge_labels[j];
        if (solutionGraph.cellOfNode(v0) == solutionGraph.cellOfNode(v1) && edge_label == 0)
        {
            auto cell = solutionGraph.cellOfNode(v0);

            auto lineProperty = graphics.defineLineProperty(
                true,
                colors[cell][0], colors[cell][1], colors[cell][2]
            );

            graphics.defineLine(v0, v1, lineProperty);
        }
        else if(edge_label == 0)
        {
            auto lineProperty = graphics.defineLineProperty(
                true,
                0, 0, 0
            );

            graphics.defineLine(v0, v1, lineProperty);
        }
    }

    graphics.center();
    graphics.normalize(0, 1);
    graphics.normalize(2);
}

} // namespace lineage

#endif
