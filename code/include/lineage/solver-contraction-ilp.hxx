#pragma once
#ifndef LINEAGE_SOLVER_CONTRACTION_ILP_HXX
#define LINEAGE_SOLVER_CONTRACTION_ILP_HXX

#include "solver-ilp.hxx"

namespace lineage {

template<typename ILP>
Solution
solveWithGraphContraction(
    Problem const& problem,
    double terminationCost = .0,
    double birthCost = .0
) {
    class DynamicGraph
    {
    public:
        DynamicGraph(size_t n) :
            vertices_(n)
        {}

        bool edgeExists(size_t a, size_t b) const
        {
            return !vertices_[a].empty() && vertices_[a].find(b) != vertices_[a].end();
        }

        std::map<size_t, double> const& getAdjacentVertices(size_t v) const
        {
            return vertices_[v];
        }

        void removeVertex(size_t v)
        {
            for (auto& p : vertices_[v])
                vertices_[p.first].erase(v);

            vertices_[v].clear();
        }

        void updateEdgeWeight(size_t a, size_t b, double w)
        {
            vertices_[a][b] += w;
            vertices_[b][a] += w;
        }

        bool vertexExists(size_t v) const
        {
            return !vertices_[v].empty();
        }

    private:
        std::vector<std::map<size_t, double>> vertices_;
    };

    andres::graph::Graph<> graph(problem.nodes.size());
    std::vector<double> edge_weights;

    Problem per_frame_problem;
    per_frame_problem.nodes = problem.nodes;
    per_frame_problem.node_offsets = problem.node_offsets;

    for (size_t i = 0; i < problem.edges.size(); ++i)
    {
        auto& e = problem.edges[i];

        if (e.t0 == e.t1)
        {
            graph.insertEdge(e.v0 + problem.node_offsets[e.t0], e.v1 + problem.node_offsets[e.t1]);
            per_frame_problem.edges.push_back(e);
        }
    }

    std::cout << "Graph before contraction:\n";
    std::cout << "..# of vertices: " << problem.nodes.size() << std::endl;
    std::cout << "..# of edges: " << problem.edges.size() << std::endl;

    auto solution = solve<andres::ilp::Gurobi<>>(per_frame_problem);

    DynamicGraph full_graph_cp(graph.numberOfVertices());

    for (size_t i = 0; i < problem.edges.size(); ++i)
    {
        auto& e = problem.edges[i];
        full_graph_cp.updateEdgeWeight(e.v0 + problem.node_offsets[e.t0], e.v1 + problem.node_offsets[e.t1], e.weight);
    }

    Problem new_problem;
    new_problem.node_offsets.push_back(0);

    size_t new_vertex_counter = 0;

    std::map<int, int> full_to_reduced_map;
    std::vector<std::vector<size_t>> partition;

    std::vector<bool> visited(graph.numberOfVertices());
    for (size_t root = 0; root < graph.numberOfVertices(); ++root)
        if (!visited[root])
        {
            std::stack<size_t> S;

            visited[root] = 1;
            S.push(root);

            if (problem.nodes[root].t == new_problem.node_offsets.size())
                new_problem.node_offsets.push_back(new_vertex_counter);

            new_problem.nodes.push_back(problem.nodes[root]);
            new_problem.nodes.back().id = new_vertex_counter - new_problem.node_offsets.back();

            full_to_reduced_map[root] = new_vertex_counter;

            partition.push_back(std::vector<size_t>());
            partition.back().push_back(root);

            while (!S.empty())
            {
                auto v = S.top();
                S.pop();

                for (auto it = graph.adjacenciesFromVertexBegin(v); it != graph.adjacenciesFromVertexEnd(v); ++it)
                {
                    auto w = it->vertex();

                    if (!visited[w] && !solution.edge_labels[it->edge()])
                    {
                        S.push(w);
                        visited[w] = 1;

                        for (auto& p : full_graph_cp.getAdjacentVertices(w))
                        {
                            if (p.first == root)
                                continue;

                            full_graph_cp.updateEdgeWeight(root, p.first, p.second);
                        }

                        full_graph_cp.removeVertex(w);
                        partition.back().push_back(w);
                    }
                }
            }

            ++new_vertex_counter;
        }

    new_problem.node_offsets.push_back(new_vertex_counter);

    andres::graph::Graph<> new_graph(new_vertex_counter);

    for (size_t v = 0; v < graph.numberOfVertices(); ++v)
        if (full_graph_cp.vertexExists(v))
            for (auto& p : full_graph_cp.getAdjacentVertices(v))
                if (!new_graph.findEdge(full_to_reduced_map[v], full_to_reduced_map[p.first]).first)
                {
                    auto& vv = new_problem.nodes[full_to_reduced_map[v]];
                    auto& ww = new_problem.nodes[full_to_reduced_map[p.first]];

                    Edge e;
                    e.t0 = vv.t;
                    e.v0 = vv.id;

                    e.t1 = ww.t;
                    e.v1 = ww.id;

                    e.weight = p.second;

                    new_problem.edges.push_back(e);

                    new_graph.insertEdge(full_to_reduced_map[v], full_to_reduced_map[p.first]);
                }

    auto sort_f = [] (Edge const& e1, Edge const& e2)
    {
        auto p1 = std::make_pair(e1.t0, e1.t1);
        auto p2 = std::make_pair(e2.t0, e2.t1);

        return p1 < p2 || (p1 == p2 && std::make_pair(e1.v0, e1.v1) < std::make_pair(e2.v0, e2.v1));
    };

    std::sort(new_problem.edges.begin(), new_problem.edges.end(), sort_f);

    std::ofstream fo("reduced-vertices.txt");

    for (auto& v : new_problem.nodes)
        fo << v.t << " " << v.id << " " << v.cx << " " << v.cy << std::endl;

    fo.close();

    fo.open("reduced-edges.txt");

    for (auto e : new_problem.edges)
        fo << e.t0 << " " << e.v0 << " " << e.t1 << " " << e.v1 << " " << e.weight << std::endl;

    fo.close();

    std::cout << "Graph after contraction:\n";
    std::cout << "..# of vertices: " << new_graph.numberOfVertices() << std::endl;
    std::cout << "..# of edges: " << new_graph.numberOfEdges() << std::endl;

    solution = solve<andres::ilp::Gurobi<>>(new_problem, terminationCost, birthCost);

    Solution full_solution;
    full_solution.node_labels.resize(problem.nodes.size());

    for (size_t i = 0; i < solution.node_labels.size(); ++i)
    {
        auto label = solution.node_labels[i];

        for (auto v : partition[i])
            full_solution.node_labels[v] = label;
    }

    full_solution.edge_labels.resize(problem.edges.size());

    for (size_t i = 0; i < problem.edges.size(); ++i)
    {
        auto& e = problem.edges[i];

        auto v0 = problem.node_offsets[e.t0] + e.v0;
        auto v1 = problem.node_offsets[e.t1] + e.v1;

        full_solution.edge_labels[i] = full_solution.node_labels[v0] == full_solution.node_labels[v1] ? 0 : 1;
    }

    return full_solution;
}

} // namespace lineage

#endif
