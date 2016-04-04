#pragma once
#ifndef LINAGE_VALIDATION_HXX
#define LINAGE_VALIDATION_HXX

#include <andres/graph/components.hxx>
#include <andres/graph/shortest-paths.hxx>

#include "problem-graph.hxx"
#include "solution.hxx"

namespace lineage {

inline
void validate(ProblemGraph const& problemGraph, Solution const& solution)
{
    struct SubgraphWithoutCut
    {
        SubgraphWithoutCut(lineage::Solution::EdgeLabels const& edgeLabels) :
            edgeLabels_(edgeLabels)
        {}
        
        bool vertex(size_t i) const
        { 
            return true;
        }
        
        bool edge(size_t i) const
        {
            return edgeLabels_[i] == 0;
        }

        lineage::Solution::EdgeLabels const& edgeLabels_;
    };

    typedef ProblemGraph::SubgraphWithoutCutAndInterFrameEdges<lineage::Solution::EdgeLabels> SubgraphCutPerFrame;
    typedef ProblemGraph::SubgraphOfTwoFramesWithoutCut<lineage::Solution::EdgeLabels> SubgraphCutTwoFrames;
    typedef andres::graph::ComponentsBySearch<ProblemGraph::Graph> ComponentLabeling;

    Problem const& problem = problemGraph.problem();

    ComponentLabeling components;
    SubgraphWithoutCut subgraphCut(solution.edge_labels);
    components.build(problemGraph.graph(), subgraphCut);

    ComponentLabeling componentsPerFrame;
    SubgraphCutPerFrame subgraphCutPerFrame(problem, solution.edge_labels);
    componentsPerFrame.build(problemGraph.graph(), subgraphCutPerFrame);

    // print numbers of nodes and edges
    std::cerr << "overall: " << problem.nodes.size()
        << " nodes, " << problem.edges.size()
        << " edges" << std::endl;

    for (size_t frame = 0; frame < problemGraph.numberOfFrames(); ++frame)
    {
        std::cerr << "frame " << frame
            << ": " << problemGraph.numberOfNodesInFrame(frame)
            << " nodes, " << problemGraph.numberOfEdgesInFrame(frame)
            << " edges" << std::endl;

        if (frame != 0)
            std::cerr << "   " << problemGraph.numberOfEdgesFromFrame(frame-1)
                << " inter-frame edges with frame " << frame - 1
                << std::endl;
        
        if (frame != problemGraph.numberOfFrames() - 1)
            std::cerr << "   " << problemGraph.numberOfEdgesFromFrame(frame)
                << " inter-frame edges with frame " << frame + 1
                << std::endl;

        // test whether edges labeled 1 well-define a multicut
        bool isMulticut = true;
        for (size_t j = 0; j < problemGraph.numberOfEdgesInFrame(frame); ++j)
        {
            auto e = problemGraph.edgeInFrame(frame, j);

            auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
            auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

            if ((solution.edge_labels[e] == 0) != (componentsPerFrame.labels_[v0] == componentsPerFrame.labels_[v1]))
            {
                std::cerr << "   error: the intra-frame edge " << e
                    << " between nodes " << v0
                    << " (t=" << problemGraph.problem().nodes[v0].t
                    << ", cx=" << problemGraph.problem().nodes[v0].cx
                    << ", cy=" << problemGraph.problem().nodes[v0].cy
                    << ", frame component=" << componentsPerFrame.labels_[v0]
                    << ") and " << v1
                    << " (t=" << problemGraph.problem().nodes[v1].t
                    << ", cx=" << problemGraph.problem().nodes[v1].cx
                    << ", cy=" << problemGraph.problem().nodes[v1].cy
                    << ", frame component=" << componentsPerFrame.labels_[v1]
                    << ") is part of a violated cycle constraint."
                    << std::endl;

                isMulticut = false;
            }
        }

        if (isMulticut)
            std::cerr << "   the intra-frame edges labeled 1 well-define a multicut" << std::endl;

        // report splits and nodes without descendants
        if (frame < problemGraph.numberOfFrames() - 1)
            for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(frame); ++j)
            {
                auto v0 = problemGraph.nodeInFrame(frame, j);

                std::set<size_t> descendantComponents;
                for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(v0); it != problemGraph.graph().adjacenciesFromVertexEnd(v0); ++it)
                {
                    auto e = it->edge();
                    auto v1 = it->vertex();

                    if (problem.nodes[v1].t == frame + 1 && solution.edge_labels[e] == 0)
                        descendantComponents.insert(componentsPerFrame.labels_[v1]);
                }

                if (descendantComponents.size() == 0)
                    std::cerr << "   warning: node " << v0
                        << " (t=" << problemGraph.problem().nodes[v0].t
                        << ", cx=" << problemGraph.problem().nodes[v0].cx
                        << ", cy=" << problemGraph.problem().nodes[v0].cy
                        << ", frame component=" << componentsPerFrame.labels_[v0]
                        << ") has *no* descendant in frame " << frame + 1
                        << std::endl;
                else if (descendantComponents.size() > 1)
                    std::cerr << "   warning: node " << v0
                        << " (t=" << problemGraph.problem().nodes[v0].cx
                        << ", cx=" << problemGraph.problem().nodes[v0].cx
                        << ", cy=" << problemGraph.problem().nodes[v0].cy
                        << ", frame component=" << componentsPerFrame.labels_[v0]
                        << ") has " << descendantComponents.size()
                        << " descendants in frame " << frame + 1
                        << std::endl;
            }

        // test morality and report nodes without ancestors
        bool isMoral = true;
        if (frame != 0)
        {
            for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(frame); ++j)
            {
                auto v0 = problemGraph.nodeInFrame(frame, j);

                std::set<size_t> ancestorComponents;
                for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(v0); it != problemGraph.graph().adjacenciesFromVertexEnd(v0); ++it)
                {
                    auto e = it->edge();
                    auto v1 = it->vertex();
                    
                    if (problem.nodes[v1].t == frame - 1 && solution.edge_labels[e] == 0)
                    {
                        ancestorComponents.insert(componentsPerFrame.labels_[v1]);

                        if (ancestorComponents.size() > 1)
                        {
                            isMoral = false;

                            std::cerr << "   error: morality violated at node " << v0
                                << " (t=" << problemGraph.problem().nodes[v0].t
                                << ", cx=" << problemGraph.problem().nodes[v0].cx
                                << ", cy=" << problemGraph.problem().nodes[v0].cy
                                << ", frame component=" << componentsPerFrame.labels_[v0]
                                << ")" << std::endl;
                        }
                    }
                }

                if (ancestorComponents.size() == 0)
                    std::cerr << "   warning: node " << v0
                        << " (t=" << problemGraph.problem().nodes[v0].t
                        << ", cx=" << problemGraph.problem().nodes[v0].cx
                        << ", cy=" << problemGraph.problem().nodes[v0].cy
                        << ", frame component=" << componentsPerFrame.labels_[v0]
                        << ") has no ancestor in frame " << frame - 1
                        << std::endl;
            }
            
            if (isMoral)
                std::cerr << "   labeling of incoming inter-frame edges is moral" << std::endl;
        }

        // test cycle constraints of inter-frame edges
        if (frame != 0)
        {
            SubgraphCutTwoFrames subgraph(problem, solution.edge_labels, frame - 1);

            bool cycleConstraintsSatisfied = true;
            for (size_t j = 0; j < problemGraph.numberOfEdgesFromFrame(frame - 1); ++j)
            {
                auto e = problemGraph.edgeFromFrame(frame - 1, j);

                if (solution.edge_labels[e] == 1) // cut
                {
                    auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = problemGraph.graph().vertexOfEdge(e, 1);
                    
                    std::deque<std::size_t> path;
                    if (andres::graph::spsp(problemGraph.graph(), subgraph, v0, v1, path) == true)
                    {
                        cycleConstraintsSatisfied = false;

                        std::cerr << "   error: the inter-frame edge " << e
                            << " between nodes " << v0
                            << " (t=" << problemGraph.problem().nodes[v0].t
                            << ", cx=" << problemGraph.problem().nodes[v0].cx
                            << ", cy=" << problemGraph.problem().nodes[v0].cy
                            << ", frame component=" << componentsPerFrame.labels_[v0]
                            << ") and " << v1
                            << " (t=" << problemGraph.problem().nodes[v1].t
                            << ", cx=" << problemGraph.problem().nodes[v1].cx
                            << ", cy=" << problemGraph.problem().nodes[v1].cy
                            << ", frame component=" << componentsPerFrame.labels_[v1]
                            << ") is part of a violated inter-frame cycle constraint."
                            << std::endl;
                    }
                }
            }

            if (cycleConstraintsSatisfied)
                std::cerr << "   labeling of incoming inter-frame edges satisfies cycle constraints" << std::endl;
        }
    }
}

} // namespace lineage

#endif
