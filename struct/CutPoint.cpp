struct CutPoint {
    int n, m, idx;
    std::vector<int> dfn, low, vis, cut;
    std::vector<std::vector<int>> adj;
    CutPoint(int _n, int _m) : n(_n), m(_m), dfn(_n + 1), 
    low(_n + 1), vis(_n + 1), cut(_n + 1), adj(_n + 1) {

    }

    void dfs(int x, int root) {
        vis[x] = 1;
        dfn[x] = ++idx;
        low[x] = idx;
        int child = 0;
        for (auto y : adj[x]) {
            if (!vis[y]) {
                dfs(y, root);
                low[x] = std::min(low[x], low[y]);
                if (low[y] >= dfn[x] && x != root) {
                    cut[x] = 1;
                }
                if (x == root) {
                    child++;
                }
            }
            low[x] = std::min(low[x], dfn[y]);
        }
        if (child >= 2 && x == root) {
            cut[x] = 1;
        }
    }

    std::vector<int> work() {
        std::vector<int> q;
        for (int i = 1; i <= n; i++) {
            if (!vis[i]) {
                dfs(i, i);
            }
        }
        for (int i = 1; i <= n; i++) {
            if (cut[i]) {
                q.push_back(i);
            }
        }
        return q;
    }

    void addEdge(int u, int v) {
        adj[u].push_back(v);
    }
};