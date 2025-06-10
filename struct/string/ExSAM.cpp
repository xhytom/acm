struct EXSAM
{
    const int CHAR_NUM = 30;   // 字符集个数，注意修改下方的 (-'a')
    int tot;                   // 节点总数：[0, tot)
    int n;
    vector<int> len, link;
    vector<vector<int>> nxt;
    EXSAM (int _n) : n(_n), len(_n * 2 + 5), link(_n * 2 + 5), nxt(n * 2 + 5, vector<int>(CHAR_NUM + 1, 0))
    {
        tot = 2;
        link[1] = -1;
    }
    int insertSAM(int last, int c)    // last 为父 c 为子
    {
        int cur = nxt[last][c];
        if (len[cur]) return cur;
        len[cur] = len[last] + 1;
        int p = link[last];
        while (p != -1)
        {
            if (!nxt[p][c])
                nxt[p][c] = cur;
            else
                break;
            p = link[p];
        }
        if (p == -1)
        {
            link[cur] = 1;
            return cur;
        }
        int q = nxt[p][c];
        if (len[p] + 1 == len[q])
        {
            link[cur] = q;
            return cur;
        }
        int clone = tot++;
        for (int i = 0; i < CHAR_NUM; ++i)
            nxt[clone][i] = len[nxt[q][i]] != 0 ? nxt[q][i] : 0;
        len[clone] = len[p] + 1;
        while (p != -1 && nxt[p][c] == q)
        {
            nxt[p][c] = clone;
            p = link[p];
        }
        link[clone] = link[q];
        link[cur] = clone;
        link[q] = clone;
        return cur;
    }

    int insertTrie(int cur, int c)
    {
        if (nxt[cur][c]) return nxt[cur][c];  // 已有该节点 直接返回
        return nxt[cur][c] = tot++;            // 无该节点 建立节点
    }

    void insert(const string &s)
    {
        int root = 1;
        for (auto ch : s) root = insertTrie(root, ch - 'a');
    }

    void insert(const char *s, int n)
    {
        int root = 1;
        for (int i = 0; i < n; ++i)
            root =
                insertTrie(root, s[i] - 'a');  // 一边插入一边更改所插入新节点的父节点
    }

    void build()
    {
        queue<pair<int, int>> q;
        for (int i = 0; i < 26; ++i)
            if (nxt[1][i]) q.push({i, 1});
        while (!q.empty())    // 广搜遍历
        {
            auto item = q.front();
            q.pop();
            auto last = insertSAM(item.second, item.first);
            for (int i = 0; i < 26; ++i)
                if (nxt[last][i]) q.push({i, last});
        }
    }
};