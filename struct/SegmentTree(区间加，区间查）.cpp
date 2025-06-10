struct SegmentTree {
    struct Node {
        int sum = 0, add = 0;
        Node() {};
    };
    int n;
    vector<int> a;
    vector<Node> tr;
    SegmentTree(int _n) : tr(_n << 2), a(_n << 2) {
        n = _n;
        build(1, 1, n);
    }
    SegmentTree(int _n, vector<int> _a) : tr(_n << 2), a(_n << 2) {
        a = _a;
        n = _n;
        build(1, 1, n);
    }
    void pushup(int u) {
        tr[u].sum = tr[u << 1].sum + tr[u << 1 | 1].sum;
    }
    void pushdown(int u, int l, int r) {
        Node &ls = tr[u << 1], &rs = tr[u << 1 | 1], &now = tr[u];
        int mid = l + r >> 1;
        ls.sum += (mid - l + 1) * now.add;
        rs.sum += (r - mid) * now.add;
        ls.add += now.add;
        rs.add += now.add;
        now.add = 0;
    }
    void build(int u, int l, int r) {
        if (l == r) {
            tr[u].sum = a[l];
            return;
        }
        int mid = l + r >> 1;
        build(u << 1, l, mid);
        build(u << 1 | 1, mid + 1, r);
        pushup(u);
    }
    void modify(int u, int L, int R, int l, int r, int x) {
        int mid = L + R >> 1;
        if (l <= L && R <= r) {
            tr[u].sum += (R - L + 1) * x;
            tr[u].add += x;
            return;
        }
        pushdown(u, L, R);       
        if (l <= mid) modify(u << 1, L, mid, l, r, x);
        if (r > mid)  modify(u << 1 | 1, mid + 1, R, l, r, x);
        pushup(u);          
    }
    int query(int u, int L, int R, int l, int r) {
        int mid = L + R >> 1;
        if (l <= L && R <= r) {
            return tr[u].sum;
        }
        pushdown(u, L, R);
        int ans = 0;
        if (l <= mid) ans += query(u << 1, L, mid, l, r);
        if (r > mid)  ans += query(u << 1 | 1, mid + 1, R, l, r);
        pushup(u);
        return ans;
    }
    void modify(int l, int r, int x) {
        modify(1, 1, n, l, r, x);
    }
    int query(int l, int r) {
        return query(1, 1, n, l, r);
    }
};