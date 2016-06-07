#pragma once
#ifndef LINEAGE_SOLUTION_GRAPH_HXX
#define LINEAGE_SOLUTION_GRAPH_HXX

#include <vector>
#include <set>
#include <map>
#include <stack>
#include <limits>

#include <andres/graph/digraph.hxx>

#include "problem-graph.hxx"
#include "validation.hxx"
#include "solution.hxx"

namespace lineage {

class SolutionGraph {
public:
    typedef andres::graph::Graph<> Graph;
    typedef andres::graph::Digraph<> Digraph;

    SolutionGraph(const ProblemGraph& problemGraph, const Solution& solution) :
        problemGraph_(problemGraph),
        solution_(solution),
        cellOfNode_(problemGraph.problem().nodes.size())
    {
        typedef andres::graph::ComponentsByPartition<Graph> ComponentLabeling;
        typedef ProblemGraph::SubgraphWithoutCutAndInterFrameEdges<Solution::EdgeLabels> SubgraphCutPerFrame;

        validate(problemGraph, solution);

        SubgraphCutPerFrame subgraphCutPerFrame(problem(), solution.edge_labels);
        ComponentLabeling componentsPerFrame;
        componentsPerFrame.build(problemGraph.graph(), subgraphCutPerFrame);

        // join components c0 (at time t) and c1 (at time t+1) iff
        // c1 is the unique descendant of c0
        for (size_t frame = 0; frame < problemGraph.numberOfFrames() - 1; ++frame)
        {
            std::map<size_t, std::set<size_t>> descendants;
            for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(frame); ++j)
            {
                auto v0 = problemGraph.nodeInFrame(frame, j);
                auto v0Component = componentsPerFrame.partition_.find(v0);
                
                for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(v0);
                    it != problemGraph.graph().adjacenciesFromVertexEnd(v0); ++it)
                {
                    auto v1 = it->vertex();

                    if (problem().nodes[v1].t == frame + 1 && solution.edge_labels[it->edge()] == 0)
                    {
                        auto v1Component = componentsPerFrame.partition_.find(v1);

                        descendants[v0Component].insert(v1Component);
                    }
                }
            }

            for (auto& p: descendants)
                if (p.second.size() == 1)
                {
                    auto v0Component = p.first;
                    auto v1Component = *(p.second.begin());

                    componentsPerFrame.partition_.merge(v0Component, v1Component);
                }
        }

        lineageGraph_.insertVertices(componentsPerFrame.partition_.numberOfSets());
        // from now on, member function numberOfCells works

        componentsPerFrame.partition_.elementLabeling(cellOfNode_.begin());

        nodesOfCell_.resize(numberOfCells());
        for (size_t v = 0; v < numberOfNodes(); ++v)
            nodesOfCell_[cellOfNode_[v]].push_back(v);

        for (size_t frame = 0; frame < problemGraph.numberOfFrames() - 1; ++frame)
            for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(frame); ++j)
            {
                auto v0 = problemGraph.nodeInFrame(frame, j);

                for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(v0);
                    it != problemGraph.graph().adjacenciesFromVertexEnd(v0); ++it)
                {
                    auto v1 = it->vertex();

                    if (problem().nodes[v1].t == frame + 1 && solution.edge_labels[it->edge()] == 0)
                    {
                        auto c0 = cellOfNode_[v0];
                        auto c1 = cellOfNode_[v1];

                        if (c0 != c1)
                            lineageGraph_.insertEdge(c0, c1);
                    }
                }
            }
    }

    size_t cellOfNode(size_t nodeIndex) const
    {
        return cellOfNode_[nodeIndex];
    }

    Digraph const& lineageGraph() const
    {
        return lineageGraph_;
    }

    size_t nodeOfCell(size_t c, size_t j) const
    {
        return nodesOfCell_[c][j];
    }

    size_t numberOfCells() const
    {
        return lineageGraph_.numberOfVertices();
    }

    size_t numberOfNodes() const
    {
        return cellOfNode_.size();
    }

    size_t numberOfNodesOfCell(size_t c) const
    {
        return nodesOfCell_[c].size();
    }

    Problem const& problem() const
    {
        return problemGraph_.problem();
    }

    ProblemGraph const& problemGraph() const
    {
        return problemGraph_;
    }

    void save(std::string const& fileNamePrefix = "lineage") const
    {
        {
            std::ofstream file(fileNamePrefix + "-fragment-node-labels.txt");
            for (size_t v = 0; v < numberOfNodes(); ++v)
                file << problem().nodes[v].t
                    << '\t' << problem().nodes[v].id
                    << '\t' << cellOfNode(v)
                    << std::endl;
            
            file.close();
        }

        saveSolution(fileNamePrefix + "-fragment-edge-labels.txt", solution_);

        {
            std::ofstream file(fileNamePrefix + "-cell-nodes.txt");
            for(size_t c = 0; c < numberOfCells(); ++c)
            {
                for(auto& v: nodesOfCell_[c])
                    file << "(" << problem().nodes[v].t
                        << " " << problem().nodes[v].id
                        << ") ";
                
                file << std::endl;
            }

            file.close();
        }
        {
            std::ofstream file(fileNamePrefix + "-cell-edges.txt");
            for(size_t c = 0; c < numberOfCells(); ++c)
                for(size_t j = 0; j < lineageGraph_.numberOfEdgesFromVertex(c); ++j)
                    file << c << '\t' << lineageGraph_.vertexFromVertex(c, j) << std::endl;

            file.close();
        }
    }

    void saveSVG(std::string const& fileName = "lineage") const
    {
        struct CellLine
        {
            size_t x;
            size_t tMin;
            size_t tMax;
        };

        std::vector<CellLine> cellLines(lineageGraph_.numberOfVertices());

        {
            std::stack<size_t> cells;

            // put all roots on the stack
            for (size_t v = 0; v < lineageGraph_.numberOfVertices(); ++v)
                if (lineageGraph_.numberOfEdgesToVertex(v) == 0)
                    cells.push(v);

            size_t xOffset = 0;
            while (!cells.empty())
            {
                auto cell = cells.top();
                cells.pop();

                cellLines[cell].x = xOffset;
                cellLines[cell].tMin = std::numeric_limits<size_t>::max();
                cellLines[cell].tMax = -std::numeric_limits<size_t>::max();

                for (auto& node : nodesOfCell_[cell])
                {
                    auto tNode = problem().nodes[node].t;
                    if (tNode < cellLines[cell].tMin)
                        cellLines[cell].tMin = tNode;
                    
                    if (tNode > cellLines[cell].tMax)
                        cellLines[cell].tMax = tNode;
                }

                if (lineageGraph_.numberOfEdgesFromVertex(cell) == 0) // if leaf
                    ++xOffset;
                else 
                    for(auto it = lineageGraph_.verticesFromVertexBegin(cell); it != lineageGraph_.verticesFromVertexEnd(cell); ++it)
                        cells.push(*it);
            }
        }

        std::ofstream file(fileName);

        // print header
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
            << std::endl
            << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
            << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
            << std::endl
            << "<svg version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
            << "xmlns:xlink=\"http://www.w3.org/1999/xlink\""
            // << " width=" << ???
            // << " height=" << ???
            << ">"
            << std::endl;

        size_t xScale = 10;
        size_t yScale = 5;
        for (size_t j = 0; j < cellLines.size(); ++j)
        {
            // draw line from (xOffset, tMin) to (xOffset, tMax)
            file << "<line x1=\"" << cellLines[j].x * xScale << "\""
                << " y1=\"" << cellLines[j].tMin * yScale<< "\""
                << " x2=\"" << cellLines[j].x * xScale << "\""
                << " y2=\"" << cellLines[j].tMax * yScale << "\""
                << " style=\"stroke:rgb(0, 0, 0); stroke-width:1pt;\"/>"
                << std::endl;

            for (size_t t = cellLines[j].tMin; t <= cellLines[j].tMax; ++t)
                file << "<line x1=\"" << (cellLines[j].x - 0.2) * xScale << "\""
                    << " y1=\"" << t * yScale<< "\""
                    << " x2=\"" << (cellLines[j].x + 0.2) * xScale << "\""
                    << " y2=\"" << t * yScale << "\""
                    << " style=\"stroke:rgb(0, 0, 0); stroke-width:1pt;\"/>"
                    << std::endl;
        }

        for (size_t e = 0; e < lineageGraph_.numberOfEdges(); ++e)
        {
            auto cell0 = lineageGraph_.vertexOfEdge(e, 0);
            auto cell1 = lineageGraph_.vertexOfEdge(e, 1);

            file << "<line x1=\"" << cellLines[cell0].x * xScale << "\""
                << " y1=\"" << cellLines[cell0].tMax * yScale<< "\""
                << " x2=\"" << cellLines[cell1].x * xScale << "\""
                << " y2=\"" << cellLines[cell1].tMin * yScale << "\""
                << " style=\"stroke:rgb(0, 153, 0); stroke-width:1pt;\"/>"
                << std::endl;
        }

        // print footer
        file << "</svg>\n";

        file.close();
    }

    Solution const& solution() const
    {
        return solution_;
    }

private:
    ProblemGraph const& problemGraph_;
    Solution const& solution_;

    std::vector<size_t> cellOfNode_;
    Digraph lineageGraph_;
    std::vector<std::vector<size_t>> nodesOfCell_;
};

} // namespace lineage

#endif
