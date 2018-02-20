#ifndef DTWCLUST_UNDIRECTEDGRAPH_HPP_
#define DTWCLUST_UNDIRECTEDGRAPH_HPP_

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dtwclust {

class UndirectedGraph {
public:
    UndirectedGraph(const int max_size);
    bool areNeighbors(const int i, const int j);
    void linkVertices(const int i, const int j);
    bool isConnected();
    int numEdges();

private:
    struct Vertex {
        Vertex(const int i) : id(i) {}
        int id;
        std::unordered_set<std::shared_ptr<Vertex>> neighbors;
    };

    int countEdges(const std::shared_ptr<Vertex>& vertex, const std::shared_ptr<Vertex>& caller);
    void dfs(const std::shared_ptr<Vertex>& vertex);

    std::unordered_map<int, std::shared_ptr<Vertex>> vertices_;
    std::vector<bool> visited_;
    int max_size_;
    bool connected_;
};

} // namespace dtwclust

#endif // DTWCLUST_UNDIRECTEDGRAPH_HPP_
