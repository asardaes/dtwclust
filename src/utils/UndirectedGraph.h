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
    void linkVertices(const int i, const int j, const int link_type);
    bool isConnected();

private:
    struct Vertex {
        Vertex(const int i, const int lnk = -1) : id(i), link_type(lnk) {}
        int id, link_type;
        std::unordered_set<std::shared_ptr<Vertex>> neighbors;
    };

    void dfs(const std::shared_ptr<Vertex>& vertex);

    std::unordered_map<int, std::shared_ptr<Vertex>> vertices_;
    std::vector<bool> visited_;
    int max_size_;
    bool connected_;
};

} // namespace dtwclust

#endif // DTWCLUST_UNDIRECTEDGRAPH_HPP_
