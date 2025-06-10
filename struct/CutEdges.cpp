struct CutEdges {
    int n;
    int idx = 0;
    vector<int> low, dfn, fa;
    vector<int> head, nxt, to, w;
    vector<int> b;
    int iddx = 1;
    vector<pair<int,int>> bridge;
    CutEdges(int n, int m) : low(n + 1), dfn(n + 1), fa(n + 1),
    head(n + 1), to(2 * m + 4), nxt(2 * m + 4), b(2 * m + 4), w(2 * m + 4), sz(n + 1, 1) {
        this->n = n;
    }
    void addEdge(int x, int y, int W) {
        nxt[++iddx] = head[x];
        head[x] = iddx;
        to[iddx] = y;
        w[iddx] = W;
    }
    vector<pair<int, int>> work() {
        for (int i = 1; i <= n; i++) {
            if (!dfn[i]) tarjan(i, 0);
        }
        return bridge;
    }
    void tarjan(int x, int e_in) {;
        dfn[x] = low[x] = ++idx;
        for(int i = head[x]; i; i = nxt[i]) {
            int y = to[i], W = w[i];
            if(!dfn[y]) {
                tarjan(y, i);
                if(dfn[x] < low[y]) {
                    bridge.push_back({x, y, W});
                    b[i] = b[i ^ 1] = 1;
                }
                low[x] = min(low[x], low[y]);
            } else if (i != (e_in ^ 1)) {
                low[x] = min(low[x], dfn[y]);
            }
        }
    }
  
};
CutEdges g(n, m);