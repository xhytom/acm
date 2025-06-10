struct SuffixAutomaton
{
    int tot, last, SIGMA = 130;
    vector<int> len, link, sz;
    vector<vector<int>> nxt;
    //vector<pii> order;
    int n;
    SuffixAutomaton(int _n) :n(_n), sz(2 * _n + 5), len(2 * _n + 5), link(2 * _n + 5), nxt(2 * _n + 5, vector<int>(SIGMA, 0))
    {
        len[1] = 0;
        link[1] = -1;
        nxt[1].clear();
        nxt[1].resize(SIGMA);
        tot = 2;
        last = 1;
    }
    void extend(int c)
    {
        int cur = tot++, p;
        len[cur] = len[last] + 1;
        nxt[cur].clear();
        nxt[cur].resize(SIGMA);
        for (p = last; p != -1 && !nxt[p][c]; p = link[p])
            nxt[p][c] = cur;
        if (p == -1) link[cur] = 1;
        else
        {
            int q = nxt[p][c];
            if (len[p] + 1 == len[q]) link[cur] = q;
            else
            {
                int clone = tot++;
                len[clone] = len[p] + 1;
                link[clone] = link[q];
                nxt[clone] = nxt[q];
                for (; p != -1 && nxt[p][c] == q; p = link[p])
                    nxt[p][c] = clone;
                link[q] = link[cur] = clone;
            }
        }
        last = cur;
        sz[cur] = 1;
    }
    vector<vector<int>> adj;
    void buildLinkTree()
    {
        adj.resize(tot + 1);
        for (int i = 2; i <= tot; i++ )
        {
            adj[link[i]].push_back(i);
        }
    }
};//sam的root为1