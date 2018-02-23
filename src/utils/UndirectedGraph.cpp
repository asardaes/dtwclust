#include "UndirectedGraph.h"

#include <algorithm> // std::fill

namespace dtwclust {

// constructor
UndirectedGraph::UndirectedGraph(const unsigned int max_size)
    : visited_(max_size, false)
    , max_size_(max_size)
    , complete_(false)
    , connected_(false)
{ }

// check if two indices are connected by an edge
bool UndirectedGraph::areNeighbors(const int i, const int j, const bool indirect) {
    auto ii = vertices_.find(i);
    if (ii == vertices_.end())
        return false;

    auto jj = vertices_.find(j);
    if (jj == vertices_.end())
        return false;

    // is it a direct neighbor?
    if (ii->second->neighbors.find(jj->second) != ii->second->neighbors.end())
        return true;

    // is it an indirect neighbor?
    if (!indirect) return false;
    std::fill(visited_.begin(), visited_.end(), false);
    this->dfs(ii->second);
    return visited_[j-1];
}

// connect two indices, creating them if they didn't exist
void UndirectedGraph::linkVertices(const int i, const int j) {
    if (i == j) return;
    std::shared_ptr<Vertex> i_vertex, j_vertex;

    auto ii = vertices_.find(i);
    if (ii == vertices_.end()) {
        i_vertex = std::make_shared<Vertex>(i);
        vertices_.insert({i, i_vertex});
    }
    else
        i_vertex = ii->second;

    auto jj = vertices_.find(j);
    if (jj == vertices_.end()) {
        j_vertex = std::make_shared<Vertex>(j);
        vertices_.insert({j, j_vertex});
    }
    else
        j_vertex = jj->second;

    i_vertex->neighbors.insert(j_vertex);
    j_vertex->neighbors.insert(i_vertex);
}

// is graph complete? https://en.wikipedia.org/wiki/Complete_graph
bool UndirectedGraph::isComplete() {
    if (complete_) return true;
    if (vertices_.size() < max_size_) return false;
    for (auto vertex : vertices_) {
        if (vertex.second->neighbors.size() != (max_size_ - 1))
            return false;
    }
    complete_ = true;
    return true;
}

// is graph conected? https://en.wikipedia.org/wiki/Connectivity_(graph_theory)
bool UndirectedGraph::isConnected() {
    if (connected_) return true;
    if (vertices_.size() < max_size_) return false;
    std::fill(visited_.begin(), visited_.end(), false);
    this->dfs(vertices_.begin()->second);
    // if all vertices were visited, the graph is connected
    // see https://www.quora.com/How-do-I-find-if-undirected-graph-is-connected-or-not-using-DFS
    for (const bool flag : visited_) {
        if (!flag) return false;
    }
    connected_ = true;
    return true;
}

// depth-first search
void UndirectedGraph::dfs(const std::shared_ptr<Vertex>& vertex) {
    // ids start with 1 due to R
    if (visited_[vertex->id - 1])
        return;
    visited_[vertex->id- 1] = true;
    for (auto neighbor : vertex->neighbors) {
        dfs(neighbor);
    }
}

} // namespace dtwclust
