#pragma once
#ifndef LINEAGE_SOLVER_LP_HXX
#define LINEAGE_SOLVER_LP_HXX

#include <cmath>
#include <vector>
#include <stack>
#include <sstream>
#include <fstream>
#include <iostream>

#include <andres/graph/shortest-paths.hxx>
#include <levinkov/timer.hxx>

#include "problem-graph.hxx"

namespace lineage {

template<typename LP>
inline
std::vector<double> solver_lp(ProblemGraph const& problemGraph, double costTermination = .0, double costBirth = .0, bool enforceBifurcationConstraint = false, std::string solution_name = "lp")
{
    const double tolerance = std::numeric_limits<float>::epsilon();

    // Dinic's Max FLow algorithm for undirected graph
    class DinicFlow
    {
    public:
        DinicFlow(int n) :
            dist_(n), q_(n), work_(n), g_(n)
        {}

        // Adds bidirectional edge
        void addEdge(int s, int t, int cap)
        {
            Edge a = { t, static_cast<int>(g_[t].size()), 0, cap };
            Edge b = { s, static_cast<int>(g_[s].size()), 0, cap };

            g_[s].push_back(a);
            g_[t].push_back(b);
        }

        void clear()
        {
            for (auto& v : g_)
                v.clear();
        }

        int maxFlow(int src, int dest)
        {
            src_ = src;
            dest_ = dest;

            for (auto& adj : g_)
                for (auto& e : adj)
                    e.f = 0;
            
            int result = 0;

            while (bfs())
            {
                std::fill(work_.begin(), work_.end(), 0);

                while (auto delta = dfs(src_, std::numeric_limits<int>::max()))
                    result += delta;
            }

            return result;
        }

        // return edges (as pairs of vertices) of the Min Cut
        std::set<std::pair<int, int>> getMinCut()
        {
            std::fill(work_.begin(), work_.end(), 0);

            std::stack<int> S;

            work_[src_] = 1;
            S.push(src_);

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                for (auto& e : g_[v])
                    if (e.f < e.cap && work_[e.to] == 0)
                    {
                        work_[e.to] = 1;
                        S.push(e.to);
                    }
            }

            std::set<std::pair<int, int>> ret;

            for (int i = 0; i < g_.size(); ++i)
                for (auto& e : g_[i])
                    if (work_[i] != work_[e.to])
                        ret.insert(std::make_pair(std::min(i, e.to), std::max(i, e.to)));

            return ret;
        }

    private:
        struct Edge
        {
            int to, rev;
            int f, cap;
        };

        bool bfs()
        {
            std::fill(dist_.begin(), dist_.end(), -1);

            dist_[src_] = 0;

            int qt = 0;
            q_[qt++] = src_;

            for (int qh = 0; qh < qt; qh++)
            {
                auto u = q_[qh];

                for (int j = 0; j < g_[u].size(); j++)
                {
                    Edge &e = g_[u][j];
                    auto v = e.to;
                    
                    if (dist_[v] < 0 && e.f < e.cap)
                    {
                        dist_[v] = dist_[u] + 1;
                        q_[qt++] = v;
                    }
                }
            }

            return dist_[dest_] >= 0;
        }

        int dfs(int u, int f)
        {
            if (u == dest_)
                return f;

            for (int &i = work_[u]; i < (int) g_[u].size(); i++)
            {
                Edge &e = g_[u][i];

                if (e.cap <= e.f) 
                    continue;

                int v = e.to;

                if (dist_[v] == dist_[u] + 1)
                {
                    auto df = dfs(v, std::min(f, e.cap - e.f));

                    if (df > 0)
                    {
                        e.f += df;
                        g_[v][e.rev].f -= df;

                        return df;
                    }
                }
            }

            return 0;
        }

        int src_, dest_;
        std::vector<int> dist_, q_, work_;
        std::vector<std::vector<Edge>> g_;
    };

    struct SubgraphOfOneFrame
    {
        SubgraphOfOneFrame(Problem const& problem) :
            problem_(problem)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            return problem_.edges[e].t0 == problem_.edges[e].t1;
        }

        Problem const& problem_;
    };

    struct SubgraphOfTwoFrames
    {
        SubgraphOfTwoFrames(Problem const& problem, size_t firstFrame) :
            problem_(problem), firstFrame_(firstFrame)
        {}

        bool vertex(size_t v) const
        {
            return true;
        }

        bool edge(size_t e) const
        {
            auto t = problem_.edges[e].t0;
            if (t != firstFrame_ && t != firstFrame_ + 1)
                return false;

            t = problem_.edges[e].t1;
            if (t != firstFrame_ && t != firstFrame_ + 1)
                return false;

            return true;
        }

        Problem const& problem_;
        size_t firstFrame_;

    };

    LP lp;
    std::vector<double> costs;

    std::vector<double> coefficients;
    std::vector<double> distances(problemGraph.graph().numberOfVertices());
    std::vector<double> edge_values(problemGraph.graph().numberOfEdges());
    std::vector<size_t> parents(problemGraph.graph().numberOfVertices());
    std::deque<size_t> path;
    std::vector<size_t> variables;

    std::ofstream file(solution_name + "-optimization-log.txt");

    size_t nIteration = 0;

    levinkov::Timer t;
    t.start();

    auto separateAndAddConstraints = [&] ()
    {
        double objValue = .0;
        for (size_t i = 0; i < costs.size(); ++i)
            objValue += costs[i]*lp.variableValue(i);

        file << t.get_elapsed_seconds() << " " << objValue;
        std::cout << t.get_elapsed_seconds() << " " << objValue << std::flush;

        levinkov::Timer t_separation;
        t_separation.start();

        DinicFlow flowInFrame(problemGraph.graph().numberOfVertices());

        for (size_t i = 0; i < edge_values.size(); ++i)
        {
            // Gurobi for some resaon sometimes produces super small but still negative values
            edge_values[i] = std::min(std::max(.0, lp.variableValue(i)), 1.0);

            auto v0 = problemGraph.graph().vertexOfEdge(i, 0);
            auto v1 = problemGraph.graph().vertexOfEdge(i, 1);

            if (problemGraph.frameOfNode(v0) == problemGraph.frameOfNode(v1))
                // take 1 - x_e in order to find a cut which is closest to satisfy to 1
                flowInFrame.addEdge(v0, v1, 100000.0*(1.0 - edge_values[i]));
        }

        size_t nSpaceCycle = 0;

        for (size_t t = 0; t < problemGraph.numberOfFrames(); ++t)
            for (size_t i = 0; i < problemGraph.numberOfEdgesInFrame(t); ++i)
            {
                auto e = problemGraph.edgeInFrame(t, i);

                auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

                if (problemGraph.frameOfNode(v0) == problemGraph.frameOfNode(v1))
                {
                    // find the shortest path using Dijkstra
                    double distance;
                    spsp(problemGraph.graph(), SubgraphOfOneFrame(problemGraph.problem()), v0, v1, edge_values.begin(), path, distance, distances.begin(), parents.begin());

                    // add inequality
                    if (edge_values[e] > distance + tolerance)
                    {
                        for (size_t j = 0; j < path.size() - 1; ++j)
                        {
                            coefficients[j] = 1.0;
                            variables[j] = problemGraph.graph().findEdge(path[j], path[j + 1]).second;
                        }

                        coefficients[path.size() - 1] = -1.0;
                        variables[path.size() - 1] = e;

                        lp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());

                        ++nSpaceCycle;
                    }
                }
            }

        file << " " << nSpaceCycle;
        std::cout << " " << nSpaceCycle << std::flush;

        auto nTotal = nSpaceCycle;
        
        size_t nSpacetimeCycle = 0;
        
        for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
            for (size_t j = 0; j < problemGraph.numberOfEdgesFromFrame(t); ++j)
            {
                auto e = problemGraph.edgeFromFrame(t, j);

                auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

                // find the shortest path using Dijkstra
                double distance;
                spsp(problemGraph.graph(), SubgraphOfTwoFrames(problemGraph.problem(), t), v0, v1, edge_values.begin(), path, distance, distances.begin(), parents.begin());

                if (edge_values[e] > distance + tolerance)
                {
                    for (size_t k = 0; k < path.size() - 1; ++k)
                    {
                        coefficients[k] = 1.0;
                        variables[k] = problemGraph.graph().findEdge(path[k], path[k + 1]).second;
                    }

                    coefficients[path.size() - 1] = -1.0;
                    variables[path.size() - 1] = e;

                    lp.addConstraint(variables.begin(), variables.begin() + path.size(), coefficients.begin(), .0, std::numeric_limits<double>::infinity());

                    ++nSpacetimeCycle;
                }
            }

        file << " " << nSpacetimeCycle;
        std::cout << " " << nSpacetimeCycle << std::flush;

        nTotal += nSpacetimeCycle;

        size_t nMorality = 0;
        
        // moral inequalities
        for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
            for (size_t i = 0; i < problemGraph.numberOfNodesInFrame(t); ++i)
                for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    if (i == j)
                        continue;

                    auto vt = problemGraph.nodeInFrame(t, i);
                    auto wt = problemGraph.nodeInFrame(t, j);

                    flowInFrame.maxFlow(vt, wt);

                    double S_cut = 1.0;

                    ptrdiff_t sz = 0;
                    for (auto& p : flowInFrame.getMinCut())
                    {
                        coefficients[sz] = -1.0;
                        variables[sz] = problemGraph.graph().findEdge(p.first, p.second).second;

                        S_cut -= 1.0 - edge_values[variables[sz]];

                        ++sz;
                    }

                    ptrdiff_t C = sz;

                    // artificial Dijkstra that searches for the shortest path that has time edges only in the beginning and the end
                    double min_distance = std::numeric_limits<double>::max();
                    std::deque<size_t> min_path;

                    for (auto vtt = problemGraph.graph().adjacenciesFromVertexBegin(vt); vtt != problemGraph.graph().adjacenciesFromVertexEnd(vt); ++vtt)
                        if (problemGraph.frameOfNode(vtt->vertex()) == t + 1)
                            for (auto wtt = problemGraph.graph().adjacenciesFromVertexBegin(wt); wtt != problemGraph.graph().adjacenciesFromVertexEnd(wt); ++wtt)
                                if (problemGraph.frameOfNode(wtt->vertex()) == t + 1)
                                {
                                    double distance;
                                    spsp(problemGraph.graph(), SubgraphOfOneFrame(problemGraph.problem()), vtt->vertex(), wtt->vertex(), edge_values.begin(), path, distance, distances.begin(), parents.begin());

                                    distance += edge_values[vtt->edge()];
                                    distance += edge_values[wtt->edge()];

                                    if (distance < min_distance)
                                    {
                                        min_distance = distance;
                                        min_path = std::move(path);
                                        min_path.push_front(vt);
                                        min_path.push_back(wt);
                                    }
                                }

                    if (S_cut > min_distance + tolerance)
                    {
                        for (size_t k = 0; k < min_path.size() - 1; ++k)
                        {
                            coefficients[sz] = 1.0;
                            variables[sz] = problemGraph.graph().findEdge(min_path[k], min_path[k + 1]).second;

                            ++sz;
                        }

                        lp.addConstraint(variables.begin(), variables.begin() + sz, coefficients.begin(), 1.0 - C, std::numeric_limits<double>::infinity());
    
                        ++nMorality;
                    }
                }

        file << " " << nMorality;
        std::cout << " " << nMorality << std::flush;

        nTotal += nMorality;

        if (costTermination > 0.0)
        {
            size_t nTermination = 0;

            for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
            {
                // build a flow network that contains edges E_t and E_{t,t+1} and a super sink
                DinicFlow flow(problemGraph.graph().numberOfVertices() + 1);

                for (size_t j = 0; j < problemGraph.numberOfEdgesInFrame(t); ++j)
                {
                    auto e = problemGraph.edgeInFrame(t, j);

                    auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = problemGraph.graph().vertexOfEdge(e, 1);
                    
                    flow.addEdge(v0, v1, 100000.0*(1.0 - edge_values[e]));
                }

                for (size_t j = 0; j < problemGraph.numberOfEdgesFromFrame(t); ++j)
                {
                    auto e = problemGraph.edgeFromFrame(t, j);

                    auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

                    flow.addEdge(v0, v1, 100000.0*(1.0 - edge_values[e]));
                }

                for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t + 1); ++j)
                    flow.addEdge(problemGraph.nodeInFrame(t + 1, j), problemGraph.graph().numberOfVertices(), std::numeric_limits<int>::max());

                for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    auto v = problemGraph.nodeInFrame(t, j);

                    flow.maxFlow(v, problemGraph.graph().numberOfVertices());

                    double S_cut = .0;

                    ptrdiff_t sz = 0;
                    for (auto& p : flow.getMinCut())
                    {
                        coefficients[sz] = -1.0;
                        variables[sz] = problemGraph.graph().findEdge(p.first, p.second).second;

                        S_cut += 1.0 - edge_values[variables[sz]];

                        ++sz;
                    }

                    if (1.0 - std::min(std::max(.0, lp.variableValue(v + problemGraph.graph().numberOfEdges())), 1.0) > S_cut + tolerance)
                    {
                        coefficients[sz] = 1.0;
                        variables[sz] = v + problemGraph.graph().numberOfEdges();

                        lp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());

                        ++nTermination;
                    }
                }
            }

            file << " " << nTermination;
            std::cout << " " << nTermination << std::flush;

            nTotal += nTermination;
        }

        if (costBirth > 0.0)
        {
            size_t nBirth = 0;

            auto offset = problemGraph.graph().numberOfEdges();
            if (costTermination > .0)
                offset += problemGraph.graph().numberOfVertices();

            for (size_t t = 1; t < problemGraph.numberOfFrames(); ++t)
            {
                // build a flow network that contains edges E_t and E_{t,t-1} and a super sink
                DinicFlow flow(problemGraph.graph().numberOfVertices() + 1);

                for (size_t j = 0; j < problemGraph.numberOfEdgesFromFrame(t - 1); ++j)
                {
                    auto e = problemGraph.edgeFromFrame(t - 1, j);

                    auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = problemGraph.graph().vertexOfEdge(e, 1);

                    flow.addEdge(v0, v1, 100000.0*(1.0 - edge_values[e]));
                }

                for (size_t j = 0; j < problemGraph.numberOfEdgesInFrame(t); ++j)
                {
                    auto e = problemGraph.edgeInFrame(t, j);

                    auto v0 = problemGraph.graph().vertexOfEdge(e, 0);
                    auto v1 = problemGraph.graph().vertexOfEdge(e, 1);
                    
                    flow.addEdge(v0, v1, 100000.0*(1.0 - edge_values[e]));
                }

                for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t - 1); ++j)
                    flow.addEdge(problemGraph.nodeInFrame(t - 1, j), problemGraph.graph().numberOfVertices(), std::numeric_limits<int>::max());

                for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t); ++j)
                {
                    auto v = problemGraph.nodeInFrame(t, j);

                    flow.maxFlow(v, problemGraph.graph().numberOfVertices());

                    double S_cut = .0;

                    int sz = 0;
                    for (auto& p : flow.getMinCut())
                    {
                        coefficients[sz] = -1.0;
                        variables[sz] = problemGraph.graph().findEdge(p.first, p.second).second;

                        S_cut += 1.0 - edge_values[variables[sz]];

                        ++sz;
                    }

                    if (1.0 - std::min(std::max(.0, lp.variableValue(v + offset)), 1.0) > S_cut + tolerance)
                    {
                        coefficients[sz] = 1.0;
                        variables[sz] = v + offset;

                        lp.addConstraint(variables.begin(), variables.begin() + sz + 1, coefficients.begin(), 1.0 - sz, std::numeric_limits<double>::infinity());

                        ++nBirth;
                    }
                }
            }

            file << " " << nBirth;
            std::cout << " " << nBirth << std::flush;

            nTotal += nBirth;
        }

        size_t nBifurcation = 0;

        if (enforceBifurcationConstraint)
        {
            size_t nBifurcation = 0;

            for (size_t t = 0; t < problemGraph.numberOfFrames() - 1; ++t)
                for (size_t i = 0; i < problemGraph.numberOfNodesInFrame(t); ++i)
                {
                    auto v = problemGraph.nodeInFrame(t, i);

                    for (size_t j = 0; j < problemGraph.numberOfNodesInFrame(t + 1); ++j)
                    {
                        auto w1 = problemGraph.nodeInFrame(t + 1, j);

                        double distance;
                        spsp(problemGraph.graph(), SubgraphOfTwoFrames(problemGraph.problem(), t), v, w1, edge_values.begin(), path, distance, distances.begin(), parents.begin());

                        if (path.size() == 0)
                            continue;

                        std::vector<size_t> path_v_w1(path.size() - 1);
                        for (size_t m = 0; m < path.size() - 1; ++m)
                            path_v_w1[m] = problemGraph.graph().findEdge(path[m], path[m + 1]).second;

                        for (size_t k = 0; k < problemGraph.numberOfNodesInFrame(t + 1); ++k)
                        {
                            if (j == k)
                                continue;

                            auto w2 = problemGraph.nodeInFrame(t + 1, k);

                            spsp(problemGraph.graph(), SubgraphOfTwoFrames(problemGraph.problem(), t), v, w2, edge_values.begin(), path, distance, distances.begin(), parents.begin());

                            if (path.size() == 0)
                                continue;

                            std::vector<size_t> path_v_w2(path.size() - 1);
                            for (size_t m = 0; m < path.size() - 1; ++m)
                                path_v_w2[m] = problemGraph.graph().findEdge(path[m], path[m + 1]).second;

                            for (size_t l = 0; l < problemGraph.numberOfNodesInFrame(t + 1); ++l)
                            {
                                if (l == j || l == k)
                                    continue;

                                auto w3 = problemGraph.nodeInFrame(t + 1, l);

                                spsp(problemGraph.graph(), SubgraphOfTwoFrames(problemGraph.problem(), t), v, w3, edge_values.begin(), path, distance, distances.begin(), parents.begin());

                                if (path.size() == 0)
                                    continue;

                                std::vector<size_t> P(path.size() - 1 + path_v_w1.size() + path_v_w2.size());

                                for (size_t m = 0; m < path.size() - 1; ++m)
                                    P[m] = problemGraph.graph().findEdge(path[m], path[m + 1]).second;

                                std::copy(path_v_w1.begin(), path_v_w1.end(), P.begin() + path.size() - 1);
                                std::copy(path_v_w2.begin(), path_v_w2.end(), P.begin() + path.size() - 1 + path_v_w1.size());

                                std::sort(P.begin(), P.end());

                                P.erase(std::unique(P.begin(), P.end()), P.end());


                                std::vector<size_t> C;

                                flowInFrame.maxFlow(w1, w2);

                                for (auto& p : flowInFrame.getMinCut())
                                    C.push_back(problemGraph.graph().findEdge(p.first, p.second).second);

                                flowInFrame.maxFlow(w1, w3);

                                for (auto& p : flowInFrame.getMinCut())
                                    C.push_back(problemGraph.graph().findEdge(p.first, p.second).second);

                                flowInFrame.maxFlow(w2, w3);

                                for (auto& p : flowInFrame.getMinCut())
                                    C.push_back(problemGraph.graph().findEdge(p.first, p.second).second);

                                std::sort(C.begin(), C.end());

                                C.erase(std::unique(C.begin(), C.end()), C.end());


                                ptrdiff_t sz = 0;

                                double S_cut = 1.0;
                                for (auto e : C)
                                {
                                    S_cut -= 1.0 - edge_values[e];

                                    coefficients[sz] = -1.0;
                                    variables[sz] = e;

                                    ++sz;
                                }

                                auto C_cut = sz;

                                distance = .0;
                                for (auto e : P)
                                {
                                    distance += edge_values[e];

                                    coefficients[sz] = 1.0;
                                    variables[sz] = e;

                                    ++sz;
                                }

                                if (S_cut > distance + tolerance)
                                {
                                    lp.addConstraint(variables.begin(), variables.begin() + sz, coefficients.begin(), 1.0 - C_cut, std::numeric_limits<double>::infinity());

                                    ++nBifurcation;
                                }
                            }
                        }
                    }
                }

            file << " " << nBifurcation;
            std::cout << " " << nBifurcation << std::flush;

            nTotal += nBifurcation;
        }

        t_separation.stop();
        t.stop();

        file << " " << t_separation.get_elapsed_seconds() << std::endl;
        std::cout << " " << t_separation.get_elapsed_seconds() << std::endl;

        std::ofstream fl(solution_name + "-variables-values-" + std::to_string(nIteration++) + ".txt");
        for (size_t i = 0; i < costs.size(); ++i)
            fl << lp.variableValue(i) << std::endl;
        fl.close();

        t.start();

        return nTotal;
    };

    for (auto e : problemGraph.problem().edges)
        costs.push_back(e.weight);

    if (costTermination > 0.0)
        costs.insert(costs.end(), problemGraph.graph().numberOfVertices(), costTermination);

    if (costBirth > 0.0)
        costs.insert(costs.end(), problemGraph.graph().numberOfVertices(), costBirth);

    coefficients.resize(costs.size());
    variables.resize(costs.size());

    lp.initModel(costs.size(), costs.data());

    for (;;)
    {
        lp.optimize();

        if (separateAndAddConstraints() == 0)
            break;
    }

    std::vector<double> solution(problemGraph.graph().numberOfEdges());
    for (size_t edge = 0; edge < problemGraph.graph().numberOfEdges(); ++edge)
        solution[edge] = lp.variableValue(edge);
}

}

#endif // of LINEAGE_SOLVER_LP_HXX