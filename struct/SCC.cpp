struct SCC {
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    int cur, cnt;
    
    SCC() {}
    SCC(int n) {
        init(n);
    }
    
    void init(int n) {
        this->n = n;
        adj.assign(n + 1, {});
        dfn.assign(n + 1, 0);
        low.resize(n + 1);
        bel.assign(n + 1, 0);
        stk.clear();
        cur = cnt = 0;
    }
    
    void addEdge(int u, int v) {
        adj[u].push_back(v);
    }
    
    void dfs(int x) {
        dfn[x] = low[x] = ++cur;
        stk.push_back(x);
        
        for (auto y : adj[x]) {
            if (dfn[y] == 0) {
                dfs(y);
                low[x] = std::min(low[x], low[y]);
            } else if (bel[y] == 0) {
                low[x] = std::min(low[x], dfn[y]);
            }
        }
        
        if (dfn[x] == low[x]) {
            int y;
            ++cnt;
            do {
                y = stk.back();
                bel[y] = cnt;
                stk.pop_back();
            } while (y != x);
        }
    }
    
    std::vector<int> work() {
        for (int i = 1; i <= n; i++) {
            if (dfn[i] == 0) {
                dfs(i);
            }
        }
        return bel;
    }
};