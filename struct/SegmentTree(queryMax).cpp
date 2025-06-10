struct SegmentTreeNode
{
    int L, R;
    int mx, add;
};

struct SegmentTree
{
    int n;
    vector<SegmentTreeNode> tr;
    SegmentTree(int _n) : tr(4 * _n), n(_n){ _build(1, 1, n);};
    void _build(int u, int l, int r)
    {
        tr[u].L = l;
        tr[u].R = r;
        int mid = l + r >> 1;
        if(l == r)
        {
            tr[u].mx = 0;
            return;
        }
        _build(u << 1, l, mid);
        _build(u << 1 | 1, mid + 1, r);
    }
    void build(int n)
    {
        _build(1, 1, n);
    }
    void _pushup(int u)
    {
        tr[u].mx = max(tr[u << 1].mx, tr[u << 1 | 1].mx);
    }
    void _pushdown(int u)
    {
        SegmentTreeNode &ls = tr[u << 1], &rs = tr[u << 1 | 1], &now = tr[u];
        ls.mx += now.add;
        ls.add += now.add;
        rs.mx += now.add;
        rs.add += now.add;
        now.add = 0;
    }
    void _modify(int u, int l, int r, int x)
    {
        int L = tr[u].L, R = tr[u].R, mid = L + R >> 1;
        if(l <= L && R <= r)
        {
            tr[u].add += x;
            tr[u].mx += x;
            return;
        }
        _pushdown(u);
        if(l <= mid)
        {
            _modify(u << 1, l, r, x);
        }
        if(r > mid)
        {
            _modify(u << 1 | 1, l, r, x);
        }
        _pushup(u);
    } 
    void modify(int l, int r, int x)
    {
        _modify(1, l, r, x);
    } 
    int _querymax(int u, int l, int r)
    {
        int L = tr[u].L, R = tr[u].R, mid = L + R >> 1;
        if(l <= L && R <= r)
        {
            return tr[u].mx;
        }
        int ans = -1e18;
        _pushdown(u);
        if(l <= mid)
        {
            ans = max(ans, _querymax(u << 1, l, r));
        }
        if(r > mid)
        {
            ans = max(ans, _querymax(u << 1 | 1, l, r));
        }
        _pushup(u);
        return ans;
    }
    int querymax(int l, int r)
    {
        return _querymax(1, l, r);
    }   
};
