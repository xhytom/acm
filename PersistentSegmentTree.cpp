struct Info {
    int val = 0;
};

struct PersistentSegmentTree {
    vector<Info> tr;
    vector<Info> a;
    vector<int> ls, rs;
    int n, idx = 1;
    PersistentSegmentTree(int _n, vector<Info> a) {
        this->n = _n;
        this->a = a;
        ls.resize(_n << 5);
        rs.resize(_n << 5);
        tr.resize(_n << 5);
        build(1, 1, _n);
    }
    void build(int u, int L, int R) {
        // test(u, L, R);
        if (L == R) {
            tr[u] = a[L];
            return;
        }
        int mid = L + R >> 1;
        if (!ls[u]) {
            ls[u] = ++idx;
        }
        if (!rs[u]) {
            rs[u] = ++idx;
        }
        build(ls[u], L, mid);
        build(rs[u], mid + 1, R);
    }   
    int modify(int u, int L, int R, int p, int x) {
        if (L == R && p == L) {
            tr[++idx] = {x};
            return idx;
        }
        int mid = L + R >> 1;
        if (p <= mid) {
            int id = modify(ls[u], L, mid, p, x);
            ls[++idx] = id;
            rs[idx] = rs[u];
            return idx;
        } else {
            int id = modify(rs[u], mid + 1, R, p, x);
            rs[++idx] = id;
            ls[idx] = ls[u];
            return idx;
        }
    }
    Info query(int u, int L, int R, int p) {
        if (L == R && L == p) {
            return tr[u];
        }
        int mid = L + R >> 1;
        if (p <= mid) {
            return query(ls[u], L, mid, p);
        } else {
            return query(rs[u], mid + 1, R, p);
        }
    }
};