struct Info {
    ll sum = 0, len = 0;
};

struct Tag {
    ll add = 0;
};

Info operator+(const Info &a, const Info &b){
    return {a.sum + b.sum, a.len + b.len};
}

void apply(Info &x, Tag &a, Tag f){
    x.sum += x.len * f.add;
    a.add += f.add;
}

template<class Info, class Tag>
struct LazySegmentTree{
    int n;
    vector<Info> info;
    vector<Tag> tag;

    LazySegmentTree() {}

    LazySegmentTree(int n, Info _init = Info()){
        init(vector<Info>(n, _init));
    }

    LazySegmentTree(const vector<Info> &_init){
        init(_init);
    }

    void init(const vector<Info> &_init){
        n = (int)_init.size() - 1;
        info.assign((n << 2) + 1, Info());
        tag.assign((n << 2) + 1, Tag());
        function<void(int, int, int)> build = [&](int p, int l, int r){
            if (l == r){
                info[p] = _init[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m + 1, r);
            pull(p);
        };
        build(1, 1, n);
    }
    
    void pull(int p){
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    
    void apply(int p, const Tag &v){
        ::apply(info[p], tag[p], v);
    }
    
    void push(int p){
        apply(2 * p, tag[p]);
        apply(2 * p + 1, tag[p]);
        tag[p] = Tag();
    }
    
    void modify(int p, int l, int r, int x, const Info &v){
        if (l == r){
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        push(p);
        if (x <= m){
            modify(2 * p, l, m, x, v);
        } 
        else{
            modify(2 * p + 1, m + 1, r, x, v);
        }
        pull(p);
    }
    
    void modify(int p, const Info &v){
        modify(1, 1, n, p, v);
    }
    
    Info query(int p, int l, int r, int x, int y){
        if (l > y || r < x){
            return Info();
        }
        if (l >= x && r <= y){
            return info[p];
        }
        int m = (l + r) / 2;
        push(p);
        return query(2 * p, l, m, x, y) + query(2 * p + 1, m + 1, r, x, y);
    }
    
    Info query(int l, int r){
        return query(1, 1, n, l, r);
    }
    
    void modify(int p, int l, int r, int x, int y, const Tag &v){
        if (l > y || r < x){
            return;
        }
        if (l >= x && r <= y){
            apply(p, v);
            return;
        }
        int m = (l + r) / 2;
        push(p);
        modify(2 * p, l, m, x, y, v);
        modify(2 * p + 1, m + 1, r, x, y, v);
        pull(p);
    }
    
    void modify(int l, int r, const Tag &v){
        return modify(1, 1, n, l, r, v);
    }
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l > y || r < x || !pred(info[p])) {
            return -1;
        }
        if (l == r) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m + 1, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 1, n, l, r, pred);
    }
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F pred) {
        if (l > y || r < x || !pred(info[p])) {
            return -1;
        }
        if (l == r) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findLast(2 * p + 1, m + 1, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 1, n, l, r, pred);
    }
};
