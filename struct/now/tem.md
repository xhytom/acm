## 有源汇上下界最大流

```c++
void solve() {
    int n, m, S, T;
    std::cin >> n >> m >> S >> T;
    std::vector<int> carry(n + 1);
    int S_ = n + 10, T_ = S_ + 1;
    g.init(S_, T_, T_);
    std::vector<std::array<int, 4>> G(m + 1);
    std::vector<int> s(m + 1);
    int ans = 0;
    for (int i = 1; i <= m; i++) {
        auto &[a, b, c, d] = G[i];
        std::cin >> a >> b >> c >> d;
        carry[a] -= c;
        carry[b] += c;
        g.addedge(a, b, d - c, i);
    }
    
    const int inf = 1e9;
    g.addedge(T, S, inf, 0);
    
    for (int i = 1; i <= n; i ++) {
        if (carry[i] > 0) {
            g.addedge(S_, i, carry[i], i + m);
            ans += carry[i];
        } else if (carry[i] < 0) {
            g.addedge(i, T_, -carry[i], i + m);
        } 
    }
    
    int D = g.dinic();
    
    g.s = S, g.t = T;
    debug(D);
    debug(ans);
    for (int i = 0; i < g.etot; i += 2) {
        auto [v, _, id, f] = g.e[i];
        auto [u, __, id_, f_] = g.e[i ^ 1];
        test(u, v, f, f_);
    }
    
    if (D >= ans) {
        std::cout << g.dinic() << '\n';
    } else {
        std::cout << "No Solution";
        return;
    }
}
```

## SCC

```c++
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
```

