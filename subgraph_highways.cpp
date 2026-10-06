#include <stdio.h>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <limits>

using namespace std;

struct HexGrid {
    int M;
    int N;
    double L;
    double width;
    double height;
    unordered_map<int, vector<double>> centers;
};

inline int cell_index(int r, int c, int M) {
    return r * M + c;
}

inline pair<int,int> index_to_cell(int index, int M) {
    return {index / M, index % M};
}

vector<int> get_hex_neighbors(int r, int c, int M, int N) {

    // Hexagonal grid neighbors (6 neighbors)

    vector<pair<int,int>> neighbors;
    if (r % 2 == 0) {
        neighbors = {{r, c-1}, {r, c+1}, {r-1, c-1}, {r-1, c}, {r+1, c-1}, {r+1, c}};
    } else {
        neighbors = {{r, c-1}, {r, c+1}, {r-1, c}, {r-1, c+1}, {r+1, c}, {r+1, c+1}};
    }
    vector<int> result;
    for (auto const& [nr, nc] : neighbors) {
        int wrapped_r = (nr % N + N) % N;
        int wrapped_c = (nc % M + M) % M;
        result.push_back(cell_index(wrapped_r, wrapped_c, M));
    }
    return result;
}

inline double periodic_dist(const vector<double>& p1, const vector<double>& p2, double width, double height) {

    // Toroidal Distance Metric

    double dx = std::abs(p1[0] - p2[0]);
    dx = std::min(dx, width - dx);
    double dy = std::abs(p1[1] - p2[1]);
    dy = std::min(dy, height - dy);
    return std::sqrt(dx*dx + dy*dy);
}

mt19937 generator(2378926);
uniform_real_distribution<double> unif(0.0, 1.0);

HexGrid make_hex_grid(double target_L) {

    // Make the hexagonal grid. Takes a target edge length and returns a grid with that approximate edge length.

    HexGrid grid;
    grid.M = max(1, (int) round(1.0 / (target_L * sqrt(3.0))));
    grid.L = 1.0 / (grid.M * sqrt(3.0));
    grid.N = max(1, (int) (2 * round((1.0 / (1.5 * grid.L)) / 2.0)));
    grid.width = 1.0;
    grid.height = grid.N * 1.5 * grid.L;

    for (int r = 0; r < grid.N; r++) {
        for (int c = 0; c < grid.M; c++) {
            double y = r * 1.5 * grid.L;
            double x = c * sqrt(3.0) * grid.L + (r % 2 == 1 ? 0.5 * sqrt(3.0) * grid.L : 0.0);
            int index = cell_index(r, c, grid.M);
            grid.centers[index] = {x, y};
        }
    }
    return grid;
}

vector<vector<double>> random_points(int n, const HexGrid& grid) {

    // Distribute n points uniformly at random in the hexagonal grid area.

    vector<vector<double>> ps;
    ps.reserve(n);
    for (int i = 0; i < n; i++) {
        ps.push_back({unif(generator) * grid.width, unif(generator) * grid.height});
    }
    return ps;
}

double s(int n) {

    // Connectivity radius (set here to be just about (log(n)/n)^1/2, the connectivity threshold)

    return pow(((log(n) + log(log(n))) / n), 0.5);
}

double p(double r, int n) {

    // Connection function, truncated so that edges longer than 3*s(n) are not included.

    if (r > 3.0 * s(n)) {
        return 0.0;
    }
    return exp(-pow(r / s(n), 2));
}

vector<vector<int>> SRGG(const vector<vector<double>>& points, const HexGrid& grid) {

    // Construct an SRGG in the specified area, on given points

    vector<vector<int>> edges;
    int n = points.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double r = periodic_dist(points[i], points[j], grid.width, grid.height);
            if (unif(generator) < p(r, n)) {
                edges.push_back({i, j});
            }
        }
    }
    return edges;
}

int point_to_cell(const vector<double>& point, const HexGrid& grid) {

    // Assigns a point to a cell in the hexagonal grid based on the closest center, using periodic distance.

    double min_dist = numeric_limits<double>::infinity();
    int min_index = -1;
    for (auto const& [cell_index, center] : grid.centers) {
        double dist = periodic_dist(point, center, grid.width, grid.height);
        if (dist < min_dist) {
            min_dist = dist;
            min_index = cell_index;
        }
    }
    return min_index;
}

vector<vector<int>> edge_list_to_adj(int n, const vector<vector<int>>& edges) {

    // Converts an edge list to an adjacency list representation of the graph.

    vector<vector<int>> adj(n);
    for (auto const& edge : edges) {
        adj[edge[0]].push_back(edge[1]);
        adj[edge[1]].push_back(edge[0]);
    }
    return adj;
}

vector<int> bfs_parents_until_targets(int n, int source, const vector<vector<int>>& adj, const vector<int>& targets) {

    // Breadth-first search from a source, until all target nodes are reached. Returns a parent array representing the BFS tree.

    vector<int> parent(n, -1);
    vector<char> is_target(n, 0);
    int remaining = 0;
    for (int t : targets) {
        if (t >= 0 && t < n) {
            is_target[t] = 1;
            remaining++;
        }
    }
    queue<int> q;
    parent[source] = source;
    q.push(source);
    if (is_target[source]) {
        remaining--;
    }
    while (!q.empty() && remaining > 0) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (parent[v] == -1) {
                parent[v] = u;
                q.push(v);
                if (is_target[v]) {
                    remaining--;
                    if (remaining == 0) break;
                }
            }
        }
    }
    return parent;
}

vector<int> weighted_parents_until_targets(
    int source,
    const vector<vector<int>>& adj,
    const vector<vector<double>>& points,
    const HexGrid& grid,
    const vector<int>& targets,
    double alpha
) {

    // Dijkstra search from a source, using periodic distance^alpha as the edge weight.

    int n = adj.size();
    vector<double> dist(n, numeric_limits<double>::infinity());
    vector<int> parent(n, -1);
    vector<char> is_target(n, 0);
    int remaining = 0;
    for (int target : targets) {
        if (target >= 0 && target < n && !is_target[target]) {
            is_target[target] = 1;
            remaining++;
        }
    }

    using PDI = pair<double, int>;
    priority_queue<PDI, vector<PDI>, greater<PDI>> pq;
    dist[source] = 0.0;
    parent[source] = source;
    pq.emplace(0.0, source);

    while (!pq.empty() && remaining > 0) {
        auto [current_dist, u] = pq.top();
        pq.pop();
        if (current_dist > dist[u]) continue;

        if (is_target[u]) {
            is_target[u] = 0;
            remaining--;
        }

        for (int v : adj[u]) {
            double edge_weight = pow(
                periodic_dist(points[u], points[v], grid.width, grid.height),
                alpha
            );
            double new_dist = dist[u] + edge_weight;
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                parent[v] = u;
                pq.emplace(new_dist, v);
            }
        }
    }
    return parent;
}


vector<pair<int,int>> build_subtree_edges(int source, const vector<int>& targets, const vector<int>& parent) {

    // Builds the edges of the subtree connecting the source to all target nodes, based on the parent array from BFS.

    vector<pair<int,int>> edges;
    vector<char> visited(parent.size(), 0);
    for (int t : targets) {
        int v = t;
        while (v != source && v != -1 && !visited[v]) {
            int p = parent[v];
            if (p == -1) break;
            edges.emplace_back(p, v);
            visited[v] = 1;
            v = p;
        }
    }
    return edges;
}

vector<int> bfs_distances(int source, const vector<vector<int>>& adj) {

    // Breadth-first search to compute distances from the source node to all other nodes in the graph represented by adjacency lists.
    // Returns a vector of distances, where distance[i] is the distance from source to node i in the BFS tree, or INT_MAX if node i is unreachable.

    int n = adj.size();
    vector<int> dist(n, INT_MAX);
    queue<int> q;
    dist[source] = 0;
    q.push(source);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (dist[v] == INT_MAX) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    return dist;
}


vector<double> minimum_energy_path(const vector<vector<double>>& points, const HexGrid& grid, const vector<vector<int>>& adj, int source, int target, double alpha=2.0, double lambda=0.0) {

    // Uses Dijkstra's algorithm to find the minimum energy path from a source to a target
    // Energy between nodes located at Z_1 and Z_2 is defined as ||Z_1-Z_2||^alpha + lambda,
    // where lambda is a regularisation parameter.

    int n = points.size();
    vector<double> dist(n, numeric_limits<double>::infinity());
    vector<int> parent(n, -1);
    dist[source] = 0.0;
    using PDI = pair<double,int>;
    priority_queue<PDI, vector<PDI>, greater<PDI>> pq;
    pq.emplace(0.0, source);
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (int v : adj[u]) {
            double energy = pow(periodic_dist(points[u], points[v], grid.width, grid.height), alpha) + lambda;

            if (dist[u] + energy < dist[v]) {
                dist[v] = dist[u] + energy;
                parent[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    vector<double> path;
    if (dist[target] == numeric_limits<double>::infinity()) {
        return path; // No path found
    }
    for (int v = target; v != -1; v = parent[v]) {
        path.push_back(v);
    }
    reverse(path.begin(), path.end());
    return path;
}


bool is_connected(const vector<vector<int>>& adj) {

    // Check whether the graph represented by adjacency lists is connected

    int n = (int)adj.size();
    if (n == 0) return true;
    vector<char> vis(n, 0);
    queue<int> q;
    q.push(0);
    vis[0] = 1;
    int seen = 1;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (!vis[v]) {
                vis[v] = 1;
                q.push(v);
                seen++;
            }
        }
    }
    return seen == n;
}


double average_stretch(const vector<vector<int>>& adj, const vector<vector<int>>& sub_adj, int alpha = 1) {

    // Computes the average stretch of the subgraph represented by sub_adj with respect to the original graph represented by adj.

    int n = adj.size();
    double total = 0.0;
    long long count = 0;
    for (int i = 0; i < n; i++) {
        vector<int> dist_base = bfs_distances(i, adj);
        vector<int> dist_sub = bfs_distances(i, sub_adj);
        for (int j = i + 1; j < n; j++) {
            if (dist_base[j] != INT_MAX && dist_sub[j] != INT_MAX) {
                total += static_cast<double>(dist_sub[j]) / dist_base[j];
                count++;
            }
        }
    }
    return (count > 0) ? total / count : numeric_limits<double>::infinity();
}

void add_adj_edge(vector<vector<int>>& adj, int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
}

vector<pair<int,int>> adjacency_to_unique_edge_list(const vector<vector<int>>& adj) {
    vector<pair<int,int>> edges;
    int n = adj.size();
    unordered_set<long long> seen;
    for (int u = 0; u < n; u++) {
        for (int v : adj[u]) {
            if (u < v) {
                long long key = ((long long)u << 32) | (unsigned long long)v;
                if (seen.insert(key).second) {
                    edges.emplace_back(u, v);
                }
            }
        }
    }
    return edges;
}

vector<vector<int>> subgraph(
    int n,
    const HexGrid& grid,
    const vector<vector<int>>& adj_list,
    const vector<vector<double>>& points,
    double alpha = 1.0
) {
    vector<vector<int>> sub_adj(n);
    unordered_map<int, vector<int>> point_assignments;
    for (int i = 0; i < n; i++) {
        int cell = point_to_cell(points[i], grid);
        point_assignments[cell].push_back(i);
    }
    unordered_map<int, int> hub_assignments;
    for (auto const& [cell, nodes] : point_assignments) {
        if (!nodes.empty()) {
            int rnd_hub = (int) floor(nodes.size() * unif(generator));
            hub_assignments[cell] = nodes[rnd_hub];
        }
    }
    for (auto const& [cell, nodes] : point_assignments) {
        if (nodes.empty()) continue;
        int hub = hub_assignments[cell];
        vector<int> parent = weighted_parents_until_targets(
            hub, adj_list, points, grid, nodes, alpha
        );
        vector<pair<int,int>> path_edges = build_subtree_edges(hub, nodes, parent);
        for (auto const& pr : path_edges) {
            add_adj_edge(sub_adj, pr.first, pr.second);
        }
    }
    for (auto const& [cell, hub] : hub_assignments) {
        auto [r, c] = index_to_cell(cell, grid.M);
        vector<int> neighbors = get_hex_neighbors(r, c, grid.M, grid.N);
        for (int neighbor_cell : neighbors) {
            if (neighbor_cell <= cell) continue;
            auto it = hub_assignments.find(neighbor_cell);
            if (it == hub_assignments.end()) continue;
            int neighbor_hub = it->second;
            vector<int> parent = weighted_parents_until_targets(
                hub, adj_list, points, grid, {neighbor_hub}, alpha
            );
            if (parent[neighbor_hub] == -1) continue;
            vector<pair<int,int>> path_edges = build_subtree_edges(hub, {neighbor_hub}, parent);
            for (auto const& pr : path_edges) {
                add_adj_edge(sub_adj, pr.first, pr.second);
            }
        }
    }
    return sub_adj;
}

int main() {
    int n = 100;
    double alpha = 2.0;
    double L = 0.11;
    HexGrid grid = make_hex_grid(L);

    // The below code generates one instance of the compression algorithm

    cout << "Generate SRGG" << endl;
    vector<vector<double>> pts = random_points(n, grid);
    vector<vector<int>> edges = SRGG(pts, grid);

    cout << "Convert to adj list" << endl;
    vector<vector<int>> adj = edge_list_to_adj(n, edges);
    cout << "Graph size: " << edges.size() << endl;

    cout << "Creating Subgraph" << endl;
    vector<vector<int>> G_hat = subgraph(n, grid, adj, pts, alpha);

    double stretch_val = average_stretch(adj, G_hat);
    vector<pair<int,int>> sub_edges = adjacency_to_unique_edge_list(G_hat);
    double edge_density = edges.empty() ? 0.0 : static_cast<double>(sub_edges.size()) / edges.size();
    cout << "Average Stretch: " << stretch_val << " | Edge Density: " << edge_density << endl;

    ofstream file("subgraph_edges.csv");
    for (auto const& edge : sub_edges) {
        file << edge.first << "," << edge.second << endl;
    }
    file.close();

    ofstream file2("graph_edges.csv");
    for (auto const& edge : edges) {
        file2 << edge[0] << "," << edge[1] << endl;
    }
    file2.close();

    ofstream file3("points.csv");
    for (auto const& point : pts) {
        file3 << point[0] << "," << point[1] << endl;
    }
    file3.close();

    // Code below runs over a range of L values to compute average stretch and edge density, but is commented out for now.

    // int M = 200;
    // double L_min = 0.2;
    // double L_max = 2.0;
    // double step = 0.1;

    // for (double L=L_min; L < L_max; L = L + step) {
    //     double stretch_avg = 0; double density_avg = 0;

    //     for (int m=0; m<M; m++) {
    //         HexGrid grid = make_hex_grid(s(n)/L);

    //         vector<vector<double>> pts = random_points(n, grid);
            
    //         vector<vector<int>> edges = SRGG(pts, grid);
    //         vector<vector<int>> adj = edge_list_to_adj(n, edges);
            
    //         vector<vector<int>> G_hat = subgraph(n, grid, adj, pts, alpha);
    //         if (!is_connected(G_hat)) {
    //             stretch_avg = stretch_avg + 10000;
    //             continue;
    //         }

    //         double stretch_val = average_stretch(adj, G_hat);
    //         vector<pair<int,int>> sub_edges = adjacency_to_unique_edge_list(G_hat);
    //         double edge_density = edges.empty() ? 0.0 : static_cast<double>(sub_edges.size()) / edges.size();

    //         stretch_avg = stretch_avg + stretch_val;
    //         density_avg = density_avg + edge_density;
    //     }

    //     cout << "L = " << L <<" Average Stretch: " << stretch_avg/static_cast<double>(M) << " | Edge Density: " << density_avg/static_cast<double>(M) << endl;
    // }



    return 0;
}
