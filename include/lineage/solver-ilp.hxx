#pragma once
#ifndef LINEAGE_SOLVER_ILP_CALLBACK_HXX
#define LINEAGE_SOLVER_ILP_CALLBACK_HXX

#include <cmath>
#include <vector>
#include <stack>
#include <sstream>
#include <fstream>
#include <iostream>

#include <andres/graph/components.hxx>
#include <andres/graph/paths.hxx>
#include <andres/graph/shortest-paths.hxx>

#include <levinkov/timer.hxx>

#include "problem-graph.hxx"
#include "solution.hxx"


namespace lineage {

template<class ILP>
Solution solver_ilp(ProblemGraph const& problemGraph, double costTermination = .0, double costBirth = .0, bool enforceBifurcationConstraint = false, std::string solutionName = "ilp")
{
    struct Data
    {
        Data(ProblemGraph const& __problemGraph) :
            problemGraph(__problemGraph)
        { }

        ProblemGraph const& problemGraph;

        double costTermination;
        double costBirth;
        std::vector<double> costs;
        bool enforceBifurcationConstraint;
        std::string solutionName;
        levinkov::Timer timer;
    };

    class Callback: public ILP::Callback
    {
    public:
        Callback(ILP& solver, Data& data) :
            ILP::Callback(solver),
            data_(data),
            coefficients_(data.costs.size()),
            variables_(data.costs.size())
        {

        }

        void separateAndAddLazyConstraints() override
        {
            std::stringstream stream;

            data_.timer.stop();
            auto time = data_.timer.get_elapsed_seconds();
            data_.timer.start();

            auto gap = (this->objectiveBest_ - this->objectiveBound_) / (1.0 + std::fabs(this->objectiveBest_));
            std::cout << time
                << ' ' << this->objectiveBound_
                << ' ' << this->objectiveBest_
                << ' ' << gap
                << std::flush;
            stream << time
                << ' ' << this->objectiveBound_
                << ' ' << this->objectiveBest_
                << ' ' << gap;

            levinkov::Timer t_separation;
            t_separation.start();

            componentsInFrame_.build(
                data_.problemGraph.graph(),
                SubgraphWithoutCutAndInterFrameEdges(data_.problemGraph.problem(), EdgeLabels(*this))
            );

            auto nSpaceCycle = separateAndAddSpaceCycleConstraints();
            std::cout << ' ' << nSpaceCycle << std::flush;
            stream << ' ' << nSpaceCycle;

            auto nSpacetimeCycle = separateAndAddSpacetimeCycleConstraints();
            std::cout << ' ' << nSpacetimeCycle << std::flush;
            stream << ' ' << nSpacetimeCycle;

            auto nMorality = separateAndAddMoralityConstraints();
            std::cout << ' ' << nMorality << std::flush;
            stream << ' ' << nMorality;

            size_t nTermination = 0;
            if (data_.costTermination > 0.0)
            {
                nTermination = separateAndAddTerminationConstraints();
                std::cout << ' ' << nTermination << std::flush;
                stream << ' ' << nTermination;
            }

            size_t nBirth = 0;
            if (data_.costBirth > 0.0)
            {
                nBirth = separateAndAddBirthConstraints();
                std::cout << ' ' << nBirth << std::flush;
                stream << ' ' << nBirth;
            }

            size_t nBifurcation = 0;
            if (data_.enforceBifurcationConstraint)
            {
                nBifurcation = separateAndAddBifurcationConstraints();
                std::cout << ' ' << nBifurcation << std::flush;
                stream << ' ' << nBifurcation;
            }

            t_separation.stop();
            data_.timer.stop(); // not keeping time for writing log

            auto objValue = .0;
            for (size_t i = 0; i < data_.costs.size(); ++i)
                objValue += data_.costs[i]*this->label(i);

            std::cout << " " << objValue;
            stream << " " << objValue;

            stream << " " << t_separation.get_elapsed_seconds() << std::endl;
            std::cout << " " << t_separation.get_elapsed_seconds() << std::endl;

            {
                std::ofstream file(data_.solutionName + "-optimization-log.txt", std::ofstream::out | std::ofstream::app);
                file << stream.str();
                file.close();
            }

            auto n = nSpaceCycle + nSpacetimeCycle + nMorality + nTermination + nBirth + nBifurcation;
            if (n == 0)
            {
                std::ofstream file(data_.solutionName + "-fragment-edge-labels-FEASIBLE-" + std::to_string(numberOfFeasibleSolutions_) + ".txt");
                for (size_t e = 0; e < data_.problemGraph.graph().numberOfEdges(); ++e)
                    file << (this->label(e) > .5 ? 1 : 0) << std::endl;
                
                file.close();

                file.open(data_.solutionName + "-variables-values-FEASIBLE-" + std::to_string(numberOfFeasibleSolutions_) + ".txt");
                for (size_t i = 0; i < data_.costs.size(); ++i)
                    file << (this->label(i) > .5 ? 1 : 0) << std::endl;
                
                file.close();

                ++numberOfFeasibleSolutions_;
            }

            ++numberOfSeparationCalls_;

            data_.timer.start(); // resume keeping time
        }

    private:
        typedef andres::graph::ComponentsBySearch<typename ProblemGraph::Graph> ComponentsType;

        class EdgeLabels
        {
        public:
            EdgeLabels(Callback& callback) :
                callback_(callback)
            {}

            int operator[](size_t edge) const
            {
                return callback_.label(edge) > .5 ? 1 : 0;
            }

        private:
            Callback& callback_;
        };

        typedef ProblemGraph::SubgraphWithoutCutAndInterFrameEdges<EdgeLabels> SubgraphWithoutCutAndInterFrameEdges;
        typedef ProblemGraph::SubgraphOfTwoFramesWithoutCut<EdgeLabels> SubgraphOfTwoFramesWithoutCut;

        size_t separateAndAddSpaceCycleConstraints()
        {
            std::deque<size_t> path;
            std::vector<ptrdiff_t> buffer;
            size_t counter = 0;

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames(); ++t)
                for (size_t i = 0; i < data_.problemGraph.numberOfEdgesInFrame(t); ++i)
                {
                    auto e = data_.problemGraph.edgeInFrame(t, i);

                    auto v0 = data_.problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = data_.problemGraph.graph().vertexOfEdge(e, 1);
                    
                    // if an edge violates connectivity as defined by connected components
                    if (
                        this->label(e) > .5
                        && data_.problemGraph.frameOfNode(v0) == data_.problemGraph.frameOfNode(v1)
                        && componentsInFrame_.areConnected(v0, v1)
                    )
                    {
                        // find the shortest path using BFS
                        spsp(
                            data_.problemGraph.graph(),
                            SubgraphWithoutCutAndInterFrameEdges(data_.problemGraph.problem(), *this),
                            v0, v1,
                            path, buffer
                        );

                        if (!findChord(data_.problemGraph.graph(), path.begin(), path.end(), true).first) // if path is chordless
                        {
                            for (size_t j = 0; j < path.size() - 1; ++j)
                            {
                                variables_[j] = data_.problemGraph.graph().findEdge(path[j], path[j + 1]).second;
                                coefficients_[j] = 1.0;
                            }

                            variables_[path.size() - 1] = e;
                            coefficients_[path.size() - 1] = -1.0;

                            this->addLazyConstraint(variables_.begin(), variables_.begin() + path.size(), coefficients_.begin(), 0, std::numeric_limits<double>::infinity());

                            ++counter;
                        }
                    }
                }

            return counter;
        }

        size_t separateAndAddSpacetimeCycleConstraints()
        {
            std::deque<size_t> path;
            std::vector<ptrdiff_t> buffer;
            size_t counter = 0;

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames() - 1; ++t)
            {
                // do connected components labeling only for frames t and t+1
                ComponentsType components;
                components.build(
                    data_.problemGraph.graph(),
                    SubgraphOfTwoFramesWithoutCut(data_.problemGraph.problem(), EdgeLabels(*this), t)
                );

                for (size_t j = 0; j < data_.problemGraph.numberOfEdgesFromFrame(t); ++j)
                {
                    auto e = data_.problemGraph.edgeFromFrame(t, j);
                    auto v0 = data_.problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = data_.problemGraph.graph().vertexOfEdge(e, 1);

                    // if a time edge violates connectivity
                    if (this->label(e) > .5 && components.areConnected(v0, v1))
                    {
                        // find the shortest path using BFS
                        spsp(
                            data_.problemGraph.graph(),
                            SubgraphOfTwoFramesWithoutCut(data_.problemGraph.problem(), EdgeLabels(*this), t),
                            v0, v1, path, buffer
                        );

                        // skip paths that have more than one time edge
                        // size_t cnt_time_edges = 0;
                        // for (auto it1 = path.begin() + 1; it1 < path.end() - 1; ++it1)
                        // {
                        //     auto v0 = *it1;
                        //     auto v1 = *(it1 + 1);

                        //     if (data_.problemGraph.frameOfNode(v0) != data_.problemGraph.frameOfNode(v1))
                        //         ++cnt_time_edges;
                        // }

                        // if (cnt_time_edges > 1)
                        //     continue;

                        // skip paths that have chords among time edges
                        bool bad = false;
                        for (auto it1 = path.begin(); it1 != path.end() - 2 && !bad; ++it1)
                            for (auto it2 = it1 + 2; it2 != path.end(); ++it2)
                            {
                                if (it1 == path.begin() && it2 == path.end() - 1)
                                    continue;

                                if (data_.problemGraph.graph().findEdge(*it1, *it2).first &&
                                    data_.problemGraph.frameOfNode(*it1) != data_.problemGraph.frameOfNode(*it2))
                                {
                                    bad = true;
                                    break;
                                }
                            }

                        if (bad)
                            continue;

                        for (size_t k = 0; k < path.size() - 1; ++k)
                        {
                            variables_[k] = data_.problemGraph.graph().findEdge(path[k], path[k + 1]).second;
                            coefficients_[k] = 1.0;
                        }

                        variables_[path.size() - 1] = e;
                        coefficients_[path.size() - 1] = -1.0;

                        this->addLazyConstraint(variables_.begin(), variables_.begin() + path.size(), coefficients_.begin(), 0, std::numeric_limits<double>::infinity());

                        ++counter;
                    }
                }
            }

            return counter;
        }

        size_t separateAndAddMoralityConstraints()
        {
            std::deque<size_t> path;
            std::vector<ptrdiff_t> buffer;
            std::vector<char> visited(data_.problemGraph.graph().numberOfVertices());
            size_t counter = 0;

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames() - 1; ++t)
            {
                // do connected components labeling only for frames t and t+1
                ComponentsType components;
                components.build(
                    data_.problemGraph.graph(),
                    SubgraphOfTwoFramesWithoutCut(data_.problemGraph.problem(), EdgeLabels(*this), t)
                );

                for (size_t i = 0; i < data_.problemGraph.numberOfNodesInFrame(t); ++i)
                    for (size_t j = 0; j < data_.problemGraph.numberOfNodesInFrame(t); ++j)
                    {
                        if (i == j)
                            continue;

                        auto v0 = data_.problemGraph.nodeInFrame(t, i);
                        auto v1 = data_.problemGraph.nodeInFrame(t, j);

                        if (components.areConnected(v0, v1) && !componentsInFrame_.areConnected(v0, v1))
                        {
                            // find the shortest path using BFS
                            spsp(
                                data_.problemGraph.graph(),
                                SubgraphOfTwoFramesWithoutCut(data_.problemGraph.problem(), EdgeLabels(*this), t),
                                v0, v1, path, buffer
                            );

                            if (path.size() < 2)
                                continue;

                            ptrdiff_t sz = 0;

                            bool bad = false;
                            for (size_t k = 0; k < path.size() - 1; ++k)
                            {
                                coefficients_[sz] = 1.0;
                                variables_[sz] = data_.problemGraph.graph().findEdge(path[k], path[k + 1]).second;

                                auto& e = data_.problemGraph.problem().edges[variables_[sz]];
                                if (e.t0 == t && e.t1 == t)
                                {
                                    bad = true;
                                    break;
                                }

                                ++sz;
                            }

                            if (bad)
                                continue;

                            // find the cut that separates v0 and v1 in frame t using DFS
                            std::fill(visited.begin() + data_.problemGraph.nodeInFrame(t, 0), visited.begin() + data_.problemGraph.nodeInFrame(t + 1, 0), 0);

                            std::stack<size_t> stack;
                            stack.push(v0);
                            visited[v0] = 1;

                            while (!stack.empty())
                            {
                                auto v = stack.top();
                                stack.pop();

                                for (auto it2 = data_.problemGraph.graph().adjacenciesFromVertexBegin(v); it2 != data_.problemGraph.graph().adjacenciesFromVertexEnd(v); ++it2)
                                {
                                    auto w = it2->vertex();

                                    if (data_.problemGraph.frameOfNode(v) == data_.problemGraph.frameOfNode(w))
                                    {
                                        if (componentsInFrame_.labels_[v] != componentsInFrame_.labels_[w])
                                        {
                                            coefficients_[sz] = -1.0;
                                            variables_[sz] = it2->edge();
                                            
                                            ++sz;
                                        }
                                        else if (!visited[w])
                                        {
                                            visited[w] = 1;
                                            stack.push(w);
                                        }
                                    }
                                }
                            }

                            // sz - path.size() + 1 = cut capacity
                            this->addLazyConstraint(variables_.begin(), variables_.begin() + sz, coefficients_.begin(), 1 - (sz - static_cast<ptrdiff_t>(path.size()) + 1), std::numeric_limits<double>::infinity());

                            ++counter;
                        }
                    }
            }

            return counter;
        }

        size_t separateAndAddTerminationConstraints()
        {
            std::vector<char> visited(data_.problemGraph.graph().numberOfVertices());
            size_t counter = 0;

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames() - 1; ++t)
                for (size_t j = 0; j < data_.problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    auto v = data_.problemGraph.nodeInFrame(t, j);
                    auto terminationVariableIndex = v + data_.problemGraph.graph().numberOfEdges();

                    // check whether the connected component in frame t is terminated or not
                    if (visited[v] == 0 && this->label(terminationVariableIndex) < .5)
                    {
                        std::vector<size_t> vertices;

                        ptrdiff_t sz = 0;
                        std::stack<size_t> stack;
                        stack.push(v);

                        visited[v] = 1;

                        bool terminated = true;

                        // do the check and find the cut at the same time
                        while (!stack.empty())
                        {
                            auto vv = stack.top();
                            stack.pop();

                            vertices.push_back(vv);

                            for (auto it = data_.problemGraph.graph().adjacenciesFromVertexBegin(vv);
                                it != data_.problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
                            {
                                auto w = it->vertex();
                                auto e = it->edge();

                                if (data_.problemGraph.frameOfNode(w) < t)
                                    continue;

                                if (data_.problemGraph.frameOfNode(w) == t + 1 && this->label(e) < .5)
                                    terminated = false;

                                if (data_.problemGraph.frameOfNode(w) == t && visited[w] == 0 && this->label(e) < .5)
                                {
                                    visited[w] = 1;
                                    stack.push(w);
                                }

                                if (this->label(e) > .5)
                                {
                                    coefficients_[sz] = -1.0;
                                    variables_[sz] = e;
                                    
                                    ++sz;
                                }
                            }
                        }

                        if (terminated)
                            // for (auto v : vertices)
                            {
                                coefficients_[sz] = 1.0;
                                variables_[sz] = v + data_.problemGraph.graph().numberOfEdges();

                                // sz = cut capacity
                                this->addLazyConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1 - sz, std::numeric_limits<double>::infinity());

                                ++counter;
                            }
                    }
                }

            return counter;
        }

        size_t separateAndAddBirthConstraints()
        {
            auto offset = data_.problemGraph.graph().numberOfEdges();
            if (data_.costTermination > .0)
                offset += data_.problemGraph.graph().numberOfVertices();

            size_t counter = 0;
            std::vector<char> visited(data_.problemGraph.graph().numberOfVertices());
            for (size_t t = 1; t < data_.problemGraph.numberOfFrames(); ++t)
                for (size_t j = 0; j < data_.problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    auto v = data_.problemGraph.nodeInFrame(t, j);
                    auto birthVariableIndex = v + offset;

                    if (visited[v] == 0 && this->label(birthVariableIndex) < .5)
                    {
                        std::vector<size_t> vertices;
                        
                        ptrdiff_t sz = 0;
                        std::stack<size_t> stack;
                        stack.push(v);
                        
                        visited[v] = 1;
                        
                        bool born = true;
                        
                        while (!stack.empty())
                        {
                            auto vv = stack.top();
                            stack.pop();
                            
                            vertices.push_back(vv);
                            
                            for (auto it = data_.problemGraph.graph().adjacenciesFromVertexBegin(vv);
                                it != data_.problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
                            {
                                auto w = it->vertex();
                                auto e = it->edge();

                                if (data_.problemGraph.frameOfNode(w) > t)
                                    continue;

                                if (data_.problemGraph.frameOfNode(w) == t - 1 && this->label(e) < .5)
                                    born = false;

                                if (data_.problemGraph.frameOfNode(w) == t && visited[w] == 0 && this->label(e) < .5)
                                {
                                    visited[w] = 1;
                                    stack.push(w);
                                }

                                if (this->label(e) > .5)
                                {
                                    coefficients_[sz] = -1.0;
                                    variables_[sz] = e;
                                    
                                    ++sz;
                                }
                            }
                        }

                        if (born)
                            // for (auto v : vertices) 
                            {
                                variables_[sz] = v + offset;
                                coefficients_[sz] = 1.0;

                                // sz = cut capacity
                                this->addLazyConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1 - sz, std::numeric_limits<double>::infinity());

                                ++counter;
                            }
                    }
                }

            return counter;
        }

        size_t separateAndAddBifurcationConstraints()
        {
            std::vector<ptrdiff_t> parent(data_.problemGraph.graph().numberOfVertices(), -1);
            size_t counter = 0;

            // for each connected component find the cut that separates it from all other components
            std::vector<std::vector<size_t>> cuts(*max_element(componentsInFrame_.labels_.begin(), componentsInFrame_.labels_.end()) + 1);

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames(); ++t)
                for (size_t i = 0; i < data_.problemGraph.numberOfEdgesInFrame(t); ++i)
                {
                    auto e = data_.problemGraph.edgeInFrame(t, i);

                    auto v0 = data_.problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = data_.problemGraph.graph().vertexOfEdge(e, 1);
                    
                    if (this->label(e) > .5)
                    {
                        cuts[componentsInFrame_.labels_[v0]].push_back(e);
                        cuts[componentsInFrame_.labels_[v1]].push_back(e);
                    }
                }

            for (auto& c : cuts)
                sort(c.begin(), c.end());

            for (size_t t = 0; t < data_.problemGraph.numberOfFrames() - 1; ++t)
                for (size_t j = 0; j < data_.problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    auto v = data_.problemGraph.nodeInFrame(t, j);

                    // if not visited
                    if (parent[v] == -1)
                    {
                        std::map<size_t, std::pair<size_t, size_t>> touch_points;

                        std::queue<size_t> Q;
                        Q.push(v);
                        
                        parent[v] = v;
                        
                        while (!Q.empty())
                        {
                            auto vv = Q.front();
                            Q.pop();
                            
                            for (auto it = data_.problemGraph.graph().adjacenciesFromVertexBegin(vv); it != data_.problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
                            {
                                auto w = it->vertex();
                                auto e = it->edge();

                                if (data_.problemGraph.frameOfNode(w) == t + 1 && this->label(e) < .5 && touch_points.count(componentsInFrame_.labels_[w]) == 0)
                                    touch_points[componentsInFrame_.labels_[w]] = std::make_pair(vv, e);

                                if (data_.problemGraph.frameOfNode(w) == t && parent[w] == -1 && this->label(e) < .5)
                                {
                                    parent[w] = vv;
                                    Q.push(w);
                                }
                            }
                        }

                        // if there are more than 2 components in the next frame
                        if (touch_points.size() > 2)
                        {
                            std::vector<size_t> P;
                            std::vector<size_t> C;

                            for (auto& p : touch_points)
                            {
                                auto w = p.second.first;
                                while (w != v)
                                {
                                    P.push_back(data_.problemGraph.graph().findEdge(parent[w], w).second);
                                    w = parent[w];
                                }

                                P.push_back(p.second.second);

                                C.insert(C.end(), cuts[p.first].begin(), cuts[p.first].end());
                            }

                            std::sort(P.begin(), P.end());
                            std::sort(C.begin(), C.end());

                            P.erase(std::unique(P.begin(), P.end()), P.end());
                            C.erase(std::unique(C.begin(), C.end()), C.end());

                            ptrdiff_t sz = 0;
                            for (auto e : P)
                            {
                                coefficients_[sz] = 1.0;
                                variables_[sz] = e;
                                
                                ++sz;
                            }

                            for (auto e : C)
                            {
                                coefficients_[sz] = -1.0;
                                variables_[sz] = e;

                                ++sz;
                            }

                            this->addLazyConstraint(variables_.begin(), variables_.begin() + sz, coefficients_.begin(), 1 - static_cast<ptrdiff_t>(C.size()), std::numeric_limits<double>::infinity());

                            ++counter;
                        }
                    }
                }

            return counter;
        }

        std::vector<double> coefficients_;
        ComponentsType components_;
        ComponentsType componentsInFrame_;
        Data& data_;
        size_t numberOfFeasibleSolutions_ { 0 };
        size_t numberOfSeparationCalls_ { 0 };
        std::vector<size_t> variables_;
    };

    // create log file/replace existing log file with empty log file
    {
        std::ofstream file(solutionName + "-optimization-log.txt");
        file.close();
    }

    Data data(problemGraph);
    data.costBirth = costBirth;
    data.costTermination = costTermination;
    data.enforceBifurcationConstraint = enforceBifurcationConstraint;
    data.solutionName = solutionName;

    // define costs
    for (auto e : problemGraph.problem().edges)
        data.costs.push_back(e.weight);
    
    if (costTermination > 0.0)
        data.costs.insert(data.costs.end(), problemGraph.graph().numberOfVertices(), costTermination);
    
    if (costBirth > 0.0)
        data.costs.insert(data.costs.end(), problemGraph.graph().numberOfVertices(), costBirth);

    ILP ilp;

    ilp.setRelativeGap(0.0);
    ilp.setAbsoluteGap(0.0);
    ilp.addVariables(data.costs.size(), data.costs.data());

    Callback callback(ilp, data);
    ilp.setCallback(callback);

    data.timer.start();
    ilp.optimize();
    data.timer.stop();

    // print runtime, objective value, bound, numbers of violated ineqs. (0)
    {
        std::stringstream stream;
        stream << data.timer.get_elapsed_seconds()
            << " " << ilp.objective()
            << " " << ilp.bound()
            << " " << ilp.gap()
            << " 0 0 0"; // violated constraints;

        if (costTermination > .0)
            stream << " 0";

        if (costBirth > .0)
            stream << " 0";

        if (enforceBifurcationConstraint)
            stream << " 0";

        stream << " 0 0\n";

        std::cout << stream.str();

        std::ofstream file(solutionName + "-optimization-log.txt", std::ofstream::out | std::ofstream::app);
        file << stream.str();
        file.close();
    }

    Solution solution;
    solution.edge_labels.resize(problemGraph.graph().numberOfEdges());
    for (size_t edge = 0; edge < problemGraph.graph().numberOfEdges(); ++edge)
        solution.edge_labels[edge] = ilp.label(edge) > 0.5 ? 1 : 0;

    return solution;
}

// template<class ILP>
// Solution solver_ilp(ProblemGraph const& problemGraph, double costTermination = .0, double costBirth = .0, bool enforceBifurcationConstraint = false, std::string solutionName = "ilp")
// {
//     typedef andres::graph::ComponentsBySearch<typename ProblemGraph::Graph> ComponentsType;

//     class EdgeLabels
//     {
//     public:
//         EdgeLabels(ILP& ilp) :
//             ilp_(ilp)
//         {}

//         double operator[](size_t edge) const
//         {
//             return ilp_.label(edge) > .5 ? 1 : 0;
//         }

//     private:
//         ILP& ilp_;
//     };

//     typedef ProblemGraph::SubgraphWithoutCutAndInterFrameEdges<EdgeLabels> SubgraphWithoutCutAndInterFrameEdges;
//     typedef ProblemGraph::SubgraphOfTwoFramesWithoutCut<EdgeLabels> SubgraphOfTwoFramesWithoutCut;

//     ILP ilp;

//     std::vector<double> costs;
//     levinkov::Timer timer;
//     std::vector<double> coefficients_;
//     ComponentsType components_;
//     ComponentsType componentsInFrame_;
//     std::vector<size_t> variables_;
    

//     auto separateAndAddSpaceCycleConstraints = [&]()
//     {
//         std::deque<size_t> path;
//         std::vector<ptrdiff_t> buffer;
//         size_t counter = 0;

//         for (size_t t = 0; t < problemGraph.numberOfFrames(); ++t)
//             for (size_t i = 0; i < problemGraph.numberOfEdgesInFrame(t); ++i)
//             {
//                 auto e = problemGraph.edgeInFrame(t, i);

//                 auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
//                 auto v1 = problemGraph.graph().vertexOfEdge(e, 1);
                
//                 // if an edge violates connectivity as defined by connected components
//                 if (
//                     ilp.label(e) > .5
//                     && problemGraph.frameOfNode(v0) == problemGraph.frameOfNode(v1)
//                     && componentsInFrame_.areConnected(v0, v1)
//                 )
//                 {
//                     // find the shortest path using BFS
//                     spsp(
//                         problemGraph.graph(),
//                         SubgraphWithoutCutAndInterFrameEdges(problemGraph.problem(), ilp),
//                         v0, v1,
//                         path, buffer
//                     );

//                     if (!findChord(problemGraph.graph(), path.begin(), path.end(), true).first) // if path is chordless
//                     {
//                         for (size_t j = 0; j < path.size() - 1; ++j)
//                         {
//                             variables_[j] = problemGraph.graph().findEdge(path[j], path[j + 1]).second;
//                             coefficients_[j] = 1.0;
//                         }

//                         variables_[path.size() - 1] = e;
//                         coefficients_[path.size() - 1] = -1.0;

//                         ilp.addConstraint(variables_.begin(), variables_.begin() + path.size(), coefficients_.begin(), 0, std::numeric_limits<double>::infinity());

//                         ++counter;
//                     }
//                 }
//             }

//         return counter;
//     };

//     auto separateAndAddSpacetimeCycleConstraints = [&]()
//     {
//         std::deque<size_t> path;
//         std::vector<ptrdiff_t> buffer;
//         size_t counter = 0;

//         for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
//         {
//             // do connected components labeling only for frames t and t+1
//             ComponentsType components;
//             components.build(
//                 problemGraph.graph(),
//                 SubgraphOfTwoFramesWithoutCut(problemGraph.problem(), EdgeLabels(ilp), t)
//             );

//             for (size_t j = 0; j < problemGraph.numberOfEdgesFromFrame(t); ++j)
//             {
//                 auto e = problemGraph.edgeFromFrame(t, j);
//                 auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
//                 auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

//                 // if a time edge violates connectivity
//                 if (ilp.label(e) > .5 && components.areConnected(v0, v1))
//                 {
//                     // find the shortest path using BFS
//                     spsp(
//                         problemGraph.graph(),
//                         SubgraphOfTwoFramesWithoutCut(problemGraph.problem(), EdgeLabels(ilp), t),
//                         v0, v1, path, buffer
//                     );

//                     // skip paths that have more than one time edge
//                     // size_t cnt_time_edges = 0;
//                     // for (auto it1 = path.begin() + 1; it1 < path.end() - 1; ++it1)
//                     // {
//                     //     auto v0 = *it1;
//                     //     auto v1 = *(it1 + 1);

//                     //     if (problemGraph.frameOfNode(v0) != problemGraph.frameOfNode(v1))
//                     //         ++cnt_time_edges;
//                     // }

//                     // if (cnt_time_edges > 1)
//                     //     continue;

//                     // skip paths that have chords among time edges
//                     bool bad = false;
//                     for (auto it1 = path.begin(); it1 != path.end() - 2 && !bad; ++it1)
//                         for (auto it2 = it1 + 2; it2 != path.end(); ++it2)
//                         {
//                             if (it1 == path.begin() && it2 == path.end() - 1)
//                                 continue;

//                             if (problemGraph.graph().findEdge(*it1, *it2).first &&
//                                 problemGraph.frameOfNode(*it1) != problemGraph.frameOfNode(*it2))
//                             {
//                                 bad = true;
//                                 break;
//                             }
//                         }

//                     if (bad)
//                         continue;

//                     for (size_t k = 0; k < path.size() - 1; ++k)
//                     {
//                         variables_[k] = problemGraph.graph().findEdge(path[k], path[k + 1]).second;
//                         coefficients_[k] = 1.0;
//                     }

//                     variables_[path.size() - 1] = e;
//                     coefficients_[path.size() - 1] = -1.0;

//                     ilp.addConstraint(variables_.begin(), variables_.begin() + path.size(), coefficients_.begin(), 0, std::numeric_limits<double>::infinity());

//                     ++counter;
//                 }
//             }
//         }

//         return counter;
//     };

//     auto separateAndAddMoralityConstraints = [&]()
//     {
//         std::deque<size_t> path;
//         std::vector<ptrdiff_t> buffer;
//         std::vector<char> visited(problemGraph.graph().numberOfVertices());
//         size_t counter = 0;

//         for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
//         {
//             // do connected components labeling only for frames t and t+1
//             ComponentsType components;
//             components.build(
//                 problemGraph.graph(),
//                 SubgraphOfTwoFramesWithoutCut(problemGraph.problem(), EdgeLabels(ilp), t)
//             );

//             for (size_t i = 0; i < problemGraph.numberOfNodesInFrame(t); ++i)
//                 for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
//                 {
//                     if (i == j)
//                         continue;

//                     auto v0 = problemGraph.nodeInFrame(t, i);
//                     auto v1 = problemGraph.nodeInFrame(t, j);

//                     if (components.areConnected(v0, v1) && !componentsInFrame_.areConnected(v0, v1))
//                     {
//                         // find the shortest path using BFS
//                         spsp(
//                             problemGraph.graph(),
//                             SubgraphOfTwoFramesWithoutCut(problemGraph.problem(), EdgeLabels(ilp), t),
//                             v0, v1, path, buffer
//                         );

//                         if (path.size() < 2)
//                             continue;

//                         ptrdiff_t sz = 0;

//                         bool bad = false;
//                         for (size_t k = 0; k < path.size() - 1; ++k)
//                         {
//                             coefficients_[sz] = 1.0;
//                             variables_[sz] = problemGraph.graph().findEdge(path[k], path[k + 1]).second;

//                             auto& e = problemGraph.problem().edges[variables_[sz]];
//                             if (e.t0 == t && e.t1 == t)
//                             {
//                                 bad = true;
//                                 break;
//                             }

//                             ++sz;
//                         }

//                         if (bad)
//                             continue;

//                         // find the cut that separates v0 and v1 in frame t using DFS
//                         std::fill(visited.begin() + problemGraph.nodeInFrame(t, 0), visited.begin() + problemGraph.nodeInFrame(t + 1, 0), 0);

//                         std::stack<size_t> stack;
//                         stack.push(v0);
//                         visited[v0] = 1;

//                         while (!stack.empty())
//                         {
//                             auto v = stack.top();
//                             stack.pop();

//                             for (auto it2 = problemGraph.graph().adjacenciesFromVertexBegin(v); it2 != problemGraph.graph().adjacenciesFromVertexEnd(v); ++it2)
//                             {
//                                 auto w = it2->vertex();

//                                 if (problemGraph.frameOfNode(v) == problemGraph.frameOfNode(w))
//                                 {
//                                     if (componentsInFrame_.labels_[v] == componentsInFrame_.labels_[w] && !visited[w])
//                                     {
//                                         visited[w] = 1;
//                                         stack.push(w);
//                                     }

//                                     if (componentsInFrame_.labels_[v] != componentsInFrame_.labels_[w])
//                                     {
//                                         coefficients_[sz] = -1.0;
//                                         variables_[sz] = it2->edge();
                                        
//                                         ++sz;
//                                     }
//                                 }
//                             }
//                         }

//                         // sz - path.size() + 1 = cut capacity
//                         ilp.addConstraint(variables_.begin(), variables_.begin() + sz, coefficients_.begin(), 1 - (sz - static_cast<ptrdiff_t>(path.size()) + 1), std::numeric_limits<double>::infinity());

//                         ++counter;
//                     }
//                 }
//         }

//         return counter;
//     };

//     auto separateAndAddTerminationConstraints = [&]()
//     {
//         std::vector<char> visited(problemGraph.graph().numberOfVertices());
//         size_t counter = 0;

//         for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
//             for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
//             {
//                 auto v = problemGraph.nodeInFrame(t, j);
//                 auto terminationVariableIndex = v + problemGraph.graph().numberOfEdges();

//                 // check whether the connected component in frame t is terminated or not
//                 if (visited[v] == 0 && ilp.label(terminationVariableIndex)  < .5)
//                 {
//                     std::vector<size_t> vertices;

//                     ptrdiff_t sz = 0;
//                     std::stack<size_t> stack;
//                     stack.push(v);

//                     visited[v] = 1;

//                     bool terminated = true;

//                     // do the check and find the cut at the same time
//                     while (!stack.empty())
//                     {
//                         auto vv = stack.top();
//                         stack.pop();

//                         vertices.push_back(vv);

//                         for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(vv);
//                             it != problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
//                         {
//                             auto w = it->vertex();
//                             auto e = it->edge();

//                             if (problemGraph.frameOfNode(w) < t)
//                                 continue;

//                             if (problemGraph.frameOfNode(w) == t + 1 && ilp.label(e) < .5)
//                                 terminated = false;

//                             if (problemGraph.frameOfNode(w) == t && visited[w] == 0 && ilp.label(e) < .5)
//                             {
//                                 visited[w] = 1;
//                                 stack.push(w);
//                             }

//                             if (ilp.label(e) > .5)
//                             {
//                                 coefficients_[sz] = -1.0;
//                                 variables_[sz] = e;
                                
//                                 ++sz;
//                             }
//                         }
//                     }

//                     if (terminated)
//                         // for (auto v : vertices)
//                         {
//                             coefficients_[sz] = 1.0;
//                             variables_[sz] = v + problemGraph.graph().numberOfEdges();

//                             // sz = cut capacity
//                             ilp.addConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1 - sz, std::numeric_limits<double>::infinity());

//                             ++counter;
//                         }
//                 }
//             }

//         return counter;
//     };

//     auto separateAndAddBirthConstraints = [&]()
//     {
//         std::vector<char> visited(problemGraph.graph().numberOfVertices());
//         size_t counter = 0;

//         for (size_t t = 1; t < problemGraph.numberOfFrames(); ++t)
//             for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
//             {
//                 auto v = problemGraph.nodeInFrame(t, j);
//                 auto birthVariableIndex = v + problemGraph.graph().numberOfEdges() + problemGraph.graph().numberOfVertices();

//                 if (visited[v] == 0 && ilp.label(birthVariableIndex) < .5)
//                 {
//                     std::vector<size_t> vertices;
                    
//                     ptrdiff_t sz = 0;
//                     std::stack<size_t> stack;
//                     stack.push(v);
                    
//                     visited[v] = 1;
                    
//                     bool born = true;
                    
//                     while (!stack.empty())
//                     {
//                         auto vv = stack.top();
//                         stack.pop();
                        
//                         vertices.push_back(vv);
                        
//                         for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(vv);
//                             it != problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
//                         {
//                             auto w = it->vertex();
//                             auto e = it->edge();

//                             if (problemGraph.frameOfNode(w) > t)
//                                 continue;

//                             if (problemGraph.frameOfNode(w) == t - 1 && ilp.label(e) < .5)
//                                 born = false;

//                             if (problemGraph.frameOfNode(w) == t && visited[w] == 0 && ilp.label(e) < .5)
//                             {
//                                 visited[w] = 1;
//                                 stack.push(w);
//                             }

//                             if (ilp.label(e) > .5)
//                             {
//                                 coefficients_[sz] = -1.0;
//                                 variables_[sz] = e;
                                
//                                 ++sz;
//                             }
//                         }
//                     }

//                     if (born)
//                         // for (auto v : vertices) 
//                         {
//                             variables_[sz] = v + problemGraph.graph().numberOfEdges() + problemGraph.graph().numberOfVertices();
//                             coefficients_[sz] = 1.0;

//                             // sz = cut capacity
//                             ilp.addConstraint(variables_.begin(), variables_.begin() + sz + 1, coefficients_.begin(), 1 - sz, std::numeric_limits<double>::infinity());

//                             ++counter;
//                         }
//                 }
//             }

//         return counter;
//     };

//     auto separateAndAddBifurcationConstraints = [&]()
//     {
//         std::vector<ptrdiff_t> parent(problemGraph.graph().numberOfVertices(), -1);
//         size_t counter = 0;

//         // for each connected component find the cut that separates it from all other components
//         std::vector<std::vector<size_t>> cuts(*max_element(componentsInFrame_.labels_.begin(), componentsInFrame_.labels_.end()) + 1);

//         for (size_t t = 0; t < problemGraph.numberOfFrames(); ++t)
//             for (size_t i = 0; i < problemGraph.numberOfEdgesInFrame(t); ++i)
//             {
//                 auto e = problemGraph.edgeInFrame(t, i);

//                 auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
//                 auto v1 = problemGraph.graph().vertexOfEdge(e, 1);
                
//                 if (ilp.label(e) > .5)
//                 {
//                     cuts[componentsInFrame_.labels_[v0]].push_back(e);
//                     cuts[componentsInFrame_.labels_[v1]].push_back(e);
//                 }
//             }

//         for (auto& c : cuts)
//             sort(c.begin(), c.end());

//         for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
//             for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
//             {
//                 auto v = problemGraph.nodeInFrame(t, j);

//                 // if not visited
//                 if (parent[v] == -1)
//                 {
//                     std::map<size_t, std::pair<size_t, size_t>> touch_points;

//                     std::queue<size_t> Q;
//                     Q.push(v);
                    
//                     parent[v] = v;
                    
//                     while (!Q.empty())
//                     {
//                         auto vv = Q.front();
//                         Q.pop();
                        
//                         for (auto it = problemGraph.graph().adjacenciesFromVertexBegin(vv); it != problemGraph.graph().adjacenciesFromVertexEnd(vv); ++it)
//                         {
//                             auto w = it->vertex();
//                             auto e = it->edge();

//                             if (problemGraph.frameOfNode(w) == t + 1 && ilp.label(e) < .5 && touch_points.count(componentsInFrame_.labels_[w]) == 0)
//                                 touch_points[componentsInFrame_.labels_[w]] = std::make_pair(vv, e);

//                             if (problemGraph.frameOfNode(w) == t && parent[w] == -1 && ilp.label(e) < .5)
//                             {
//                                 parent[w] = vv;
//                                 Q.push(w);
//                             }
//                         }
//                     }

//                     // if there are more than 2 components in the next frame
//                     if (touch_points.size() > 2)
//                     {
//                         std::vector<size_t> P;
//                         std::vector<size_t> C;

//                         for (auto& p : touch_points)
//                         {
//                             auto w = p.second.first;
//                             while (w != v)
//                             {
//                                 P.push_back(problemGraph.graph().findEdge(parent[w], w).second);
//                                 w = parent[w];
//                             }

//                             P.push_back(p.second.second);

//                             C.insert(C.end(), cuts[p.first].begin(), cuts[p.first].end());
//                         }

//                         std::sort(P.begin(), P.end());
//                         std::sort(C.begin(), C.end());

//                         P.erase(std::unique(P.begin(), P.end()), P.end());
//                         C.erase(std::unique(C.begin(), C.end()), C.end());

//                         ptrdiff_t sz = 0;
//                         for (auto e : P)
//                         {
//                             coefficients_[sz] = 1.0;
//                             variables_[sz] = e;
                            
//                             ++sz;
//                         }

//                         for (auto e : C)
//                         {
//                             coefficients_[sz] = -1.0;
//                             variables_[sz] = e;

//                             ++sz;
//                         }

//                         ilp.addConstraint(variables_.begin(), variables_.begin() + sz, coefficients_.begin(), 1 - static_cast<ptrdiff_t>(C.size()), std::numeric_limits<double>::infinity());

//                         ++counter;
//                     }
//                 }
//             }

//         return counter;
//     };

//     auto separateAndAddLazyConstraints = [&]()
//     {
//         std::stringstream stream;

//         timer.stop();
//         auto time = timer.get_elapsed_seconds();
//         timer.start();

//         std::cout << time << std::flush;
//         stream << time;

//         levinkov::Timer t_separation;
//         t_separation.start();

//         componentsInFrame_.build(
//             problemGraph.graph(),
//             SubgraphWithoutCutAndInterFrameEdges(problemGraph.problem(), EdgeLabels(ilp))
//         );

//         auto nSpaceCycle = separateAndAddSpaceCycleConstraints();
//         std::cout << ' ' << nSpaceCycle << std::flush;
//         stream << ' ' << nSpaceCycle;

//         auto nSpacetimeCycle = separateAndAddSpacetimeCycleConstraints();
//         std::cout << ' ' << nSpacetimeCycle << std::flush;
//         stream << ' ' << nSpacetimeCycle;

//         auto nMorality = separateAndAddMoralityConstraints();
//         std::cout << ' ' << nMorality << std::flush;
//         stream << ' ' << nMorality;

//         size_t nTermination = 0;
//         if (costTermination > 0.0)
//         {
//             nTermination = separateAndAddTerminationConstraints();
//             std::cout << ' ' << nTermination << std::flush;
//             stream << ' ' << nTermination;
//         }

//         size_t nBirth = 0;
//         if (costBirth > 0.0)
//         {
//             nBirth = separateAndAddBirthConstraints();
//             std::cout << ' ' << nBirth << std::flush;
//             stream << ' ' << nBirth;
//         }

//         size_t nBifurcation = 0;
//         if (enforceBifurcationConstraint)
//         {
//             nBifurcation = separateAndAddBifurcationConstraints();
//             std::cout << ' ' << nBifurcation << std::flush;
//             stream << ' ' << nBifurcation;
//         }

//         t_separation.stop();
//         timer.stop(); // not keeping time for writing log

//         auto objValue = .0;
//         for (size_t i = 0; i < costs.size(); ++i)
//             objValue += costs[i]*ilp.label(i);

//         std::cout << " " << objValue;
//         stream << " " << objValue;

//         stream << " " << t_separation.get_elapsed_seconds() << std::endl;
//         std::cout << " " << t_separation.get_elapsed_seconds() << std::endl;

//         {
//             std::ofstream file(solutionName + "-optimization-log.txt", std::ofstream::out | std::ofstream::app);
//             file << stream.str();
//             file.close();
//         }

//         auto n = nSpaceCycle + nSpacetimeCycle + nMorality + nTermination + nBirth + nBifurcation;

//         // ++numberOfSeparationCalls_;

//         timer.start(); // resume keeping time

//         return n;
//     };

//     // create log file/replace existing log file with empty log file
//     {
//         std::ofstream file(solutionName + "-optimization-log.txt");
//         file.close();
//     }

//     // define costs
//     for (auto e : problemGraph.problem().edges)
//         costs.push_back(e.weight);

//     if (costTermination > 0.0)
//         costs.insert(costs.end(), problemGraph.graph().numberOfVertices(), costTermination);
    
//     if (costBirth > 0.0)
//         costs.insert(costs.end(), problemGraph.graph().numberOfVertices(), costBirth);


//     variables_.resize(costs.size());
//     coefficients_.resize(costs.size());

//     ilp.initModel(costs.size(), costs.data());

//     timer.start();
//     size_t n = 1;
//     while (n)
//     {
//         ilp.optimize();

//         n = separateAndAddLazyConstraints();
//     }
//     timer.stop();

//     // print runtime, objective value, bound, numbers of violated ineqs. (0)
//     // {
//     //     std::stringstream stream;
//     //     stream << timer.get_elapsed_seconds()
//     //         << " " << ilp.objective()
//     //         << " " << ilp.bound()
//     //         << " " << ilp.gap()
//     //         << " 0 0 0 0 0" // violated constraints
//     //         << std::endl;

//     //     std::cout << stream.str();

//     //     std::ofstream file(solutionName + "-optimization-log.txt", std::ofstream::out | std::ofstream::app);
//     //     file << stream.str();
//     //     file.close();
//     // }

//     // save solution
//     std::ofstream fl(solutionName + "-variables-values.txt");
//     for (int i = 0; i < costs.size(); ++i)
//         fl << ilp.label(i) << std::endl;
//     fl.close();

//     double objValue = .0;
//     // for (int i = 0; i < costs.size(); ++i)
//     for (int i = 0; i < problemGraph.graph().numberOfEdges(); ++i)
//         objValue += costs[i]*ilp.label(i);
//     std::cout << "Objective value: " << objValue << std::endl;

//     Solution solution;
//     solution.edge_labels.resize(problemGraph.graph().numberOfEdges());
//     for (size_t edge = 0; edge < problemGraph.graph().numberOfEdges(); ++edge)
//         solution.edge_labels[edge] = ilp.label(edge) > 0.5 ? 1 : 0;

//     return solution;
// }

} // namespace lineage

#endif
