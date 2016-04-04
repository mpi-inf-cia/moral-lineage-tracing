#pragma once
#ifndef LINEAGE_PROBLEM_GRAPH_HXX
#define LINEAGE_PROBLEM_GRAPH_HXX

#include <vector>

#include <andres/graph/graph.hxx>

#include "problem.hxx"

namespace lineage {

class ProblemGraph {
public:
    typedef andres::graph::Graph<> Graph;

    template<class EdgeLabels>
    struct SubgraphWithoutCutAndInterFrameEdges
    {
        SubgraphWithoutCutAndInterFrameEdges(Problem const& problem, EdgeLabels const& edgeLabels) :
            problem_(problem), edgeLabels_(edgeLabels)
        {}

        bool vertex(size_t i) const
        {
            return true;
        }

        bool edge(size_t i) const
        {
            return edgeLabels_[i] == 0 && problem_.edges[i].t0 == problem_.edges[i].t1;
        }

        Problem const& problem_;
        EdgeLabels const& edgeLabels_;
    };

    template<class EdgeLabels>
    struct SubgraphOfTwoFramesWithoutCut
    {
        SubgraphOfTwoFramesWithoutCut(Problem const& problem, EdgeLabels const& edgeLabels, size_t firstFrame) :
            problem_(problem), edgeLabels_(edgeLabels), firstFrame_(firstFrame)
        {}

        bool vertex(size_t i) const
        {
            auto t = problem_.nodes[i].t;
            
            return t == firstFrame_ || t == firstFrame_ + 1;
        }

        bool edge(size_t i) const
        {
            if(edgeLabels_[i] != 0)
                return false;

            auto t = problem_.edges[i].t0;
            if (t != firstFrame_ && t != firstFrame_ + 1)
                return false;

            t = problem_.edges[i].t1;
            if (t != firstFrame_ && t != firstFrame_ + 1)
                return false;

            return true;
        }

        Problem const& problem_;
        EdgeLabels const& edgeLabels_;
        size_t firstFrame_;
    };

    ProblemGraph(Problem const& problem) :
        problem_(problem)
    {
        for(size_t j = 0; j < problem.nodes.size(); ++j)
        {
            auto const& node = problem.nodes[j];
            if (node.t >= nodeIndicesInFrame_.size())
                nodeIndicesInFrame_.resize(node.t + 1);

            nodeIndicesInFrame_[node.t].push_back(j);
        }
        
        numberOfFrames_ = nodeIndicesInFrame_.size();

        graph_.insertVertices(problem.nodes.size());

        edgeIndicesInFrame_.resize(numberOfFrames_);
        edgeIndicesFromFrame_.resize(numberOfFrames_ - 1);

        for (size_t j = 0; j < problem.edges.size(); ++j)
        {
            auto const& edge = problem.edges[j];

            if (edge.t0 == edge.t1)
                edgeIndicesInFrame_[edge.t0].push_back(j);
            else
            {
                if (edge.t0 + 1 != edge.t1)
                    throw std::runtime_error("edge is not directed to next frame.");

                edgeIndicesFromFrame_[edge.t0].push_back(j);
            }

            graph_.insertEdge(nodeIndicesInFrame_[edge.t0][edge.v0], nodeIndicesInFrame_[edge.t1][edge.v1]);
        }
    }


    size_t edgeFromFrame(size_t t, size_t j) const
    {
        return edgeIndicesFromFrame_[t][j]; 
    }

    size_t edgeInFrame(size_t t, size_t j) const
    {
        return edgeIndicesInFrame_[t][j];
    }

    size_t frameOfNode(size_t v) const
    {
        return problem_.nodes[v].t;
    }

    Graph const& graph() const
    {
        return graph_;
    }

    size_t numberOfEdgesFromFrame(size_t t) const
    {
        return edgeIndicesFromFrame_[t].size();
    }

    size_t numberOfEdgesInFrame(size_t t) const
    {
        return edgeIndicesInFrame_[t].size();
    }

    size_t numberOfFrames() const
    {
        return numberOfFrames_;
    }

    size_t numberOfNodesInFrame(size_t t) const
    {
        return nodeIndicesInFrame_[t].size();
    }

    size_t nodeInFrame(size_t t, size_t j) const
    {
        return nodeIndicesInFrame_[t][j];
    }

    Problem const& problem() const
    {
        return problem_;
    }

private:
    Problem const& problem_;

    Graph graph_;
    std::vector<std::vector<size_t>> edgeIndicesFromFrame_;
    std::vector<std::vector<size_t>> edgeIndicesInFrame_;
    std::vector<std::vector<size_t>> nodeIndicesInFrame_;
    size_t numberOfFrames_;
};

} // namespace lineage

#endif
