struct HLD {
    int n, root;
    std::vector<std::vector<int>> adj;
    std::vector<int> dep, sz, hs, top, fa;
    HLD (int _n, int _root, std::vector<std::vector<int>> &_adj) : n(_n), root(_root), 
        dep(n + 1), sz(n + 1), hs(n + 1), top(n + 1), fa(n + 1) {
        adj = _adj;
        dfs0(root, 0);
        dfs1(root, 0, root);
    }
    void dfs0 (int x, int fx) {
        sz[x] = 1;
        fa[x] = fx;
        dep[x] = dep[fx] + 1;
        for (auto y : adj[x]) {
            if (y != fx) {
                dfs0(y, x);
                sz[x] += sz[y];
                if (sz[y] > sz[hs[x]]) {
                    hs[x] = y;
                }
            }
        }
    }
    void dfs1 (int x, int fx, int root) {
        top[x] = root;
        if (hs[x]) {
            dfs1(hs[x], x, root);
        }
        for (auto y : adj[x]) {
            if (y != fx && y != hs[x]) {
                dfs1(y, x, y);
            }
        }
    }
    int getLca (int x, int y) {
        while (top[x] != top[y]) {
            if (dep[top[x]] < dep[top[y]]) swap(x, y);
            x = fa[top[x]];
        }
        if (dep[x] < dep[y]) {
            return x;
        } else {
            return y;
        }
    }
};