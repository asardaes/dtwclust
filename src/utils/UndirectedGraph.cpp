#include "UndirectedGraph.h"

#include <algorithm> // std::fill
#include <cstddef> // std::size_t
#include <functional> // std::hash
#include <memory>
#include <unordered_set>

namespace dtwclust {

// =================================================================================================
/* Vertex
 *
 * See https://stackoverflow.com/q/27085782/5793905 and
 * https://stackoverflow.com/questions/13695640/how-to-make-a-c11-stdunordered-set-of-stdweak-ptr
 *
 * This should be fine in my case because vertices_ always has all existing shared_ptr<Vertex>,
 * and it will never go out of scope before neighbors.
 */

// struct definition
struct UndirectedGraph::Vertex {
    Vertex(const int i) : id(i) {}
    int id;
    std::unordered_set<std::weak_ptr<Vertex>, VertexHash, VertexEqual> neighbors;
};

// hash
std::size_t UndirectedGraph::VertexHash::operator()(const std::weak_ptr<Vertex>& vertex_wp) const {
    return std::hash<int>{}(vertex_wp.lock()->id);
}

// equal_to
bool UndirectedGraph::VertexEqual::operator()(const std::weak_ptr<Vertex>& left,
                                              const std::weak_ptr<Vertex>& right) const
{
    return left.lock()->id == right.lock()->id;
}

// =================================================================================================
/* UndirectedGraph */

// constructor
UndirectedGraph::UndirectedGraph(const unsigned int max_size)
    : visited_(max_size, false)
    , max_size_(max_size)
    , complete_(false)
    , connected_(false)
{ }

// check if two indices are connected directly by an edge
bool UndirectedGraph::areNeighbors(const int i, const int j) {
    auto ii = vertices_.find(i);
    if (ii == vertices_.end())
        return false;

    auto jj = vertices_.find(j);
    if (jj == vertices_.end())
        return false;

    // is it a direct neighbor?
    if (ii->second->neighbors.find(jj->second) != ii->second->neighbors.end())
        return true;
    return false;
}

// connect two indices, creating them if they didn't exist
void UndirectedGraph::linkVertices(const int i, const int j, const bool deeply) {
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

    if (deeply) {
        std::fill(visited_.begin(), visited_.end(), false);
        this->dfs(i_vertex);
        int size = visited_.size();
        for (int i_visited = 1; i_visited < size; i_visited++) {
            for (int j_visited = 0; j_visited < i_visited; j_visited++) {
                if (visited_[i_visited] && visited_[j_visited])
                    this->linkVertices(i_visited + 1, j_visited + 1, false);
            }
        }
    }
}

// is graph complete? https://en.wikipedia.org/wiki/Complete_graph
bool UndirectedGraph::isComplete() {
    if (complete_) return true;
    if (vertices_.size() < max_size_) return false;
    for (auto vertex : vertices_) {
        if (vertex.second->neighbors.size() != (max_size_ - 1)) return false;
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
    for (const bool flag : visited_) if (!flag) return false;
    connected_ = true;
    return true;
}

// depth-first search
void UndirectedGraph::dfs(const std::shared_ptr<Vertex>& vertex) {
    // ids start with 1 due to R
    if (visited_[vertex->id - 1])
        return;
    visited_[vertex->id - 1] = true;
    for (auto neighbor : vertex->neighbors)
        dfs(neighbor.lock());
}

} // namespace dtwclust
