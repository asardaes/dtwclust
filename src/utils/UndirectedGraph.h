#ifndef DTWCLUST_UNDIRECTEDGRAPH_HPP_
#define DTWCLUST_UNDIRECTEDGRAPH_HPP_

#include <cstddef> // std::size_t
#include <memory>
#include <unordered_map>
#include <vector>

namespace dtwclust {

class UndirectedGraph {
public:
    UndirectedGraph(const unsigned int max_size);
    bool areNeighbors(const int i, const int j);
    void linkVertices(const int i, const int j, const bool deeply = false);
    bool isComplete();
    bool isConnected();

    struct Vertex;

private:
    void dfs(const std::shared_ptr<Vertex>& vertex);

    std::unordered_map<int, std::shared_ptr<Vertex>> vertices_;
    std::vector<bool> visited_;
    unsigned int max_size_;
    bool complete_, connected_;
};

} // namespace dtwclust

#endif // DTWCLUST_UNDIRECTEDGRAPH_HPP_
