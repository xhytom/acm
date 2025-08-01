struct Info {
    
};
 
Info operator+(const Info &a, const Info &b){
    
}
 
template<class Info>
struct SegmentTree{
    int n;
    vector<Info> info;
 
    SegmentTree() {}
 
    SegmentTree(int n, Info _init = Info()){
        init(vector<Info>(n, _init));
    }
 
    SegmentTree(const vector<Info> &_init){
        init(_init);
    }
 
    void init(const vector<Info> &_init){
        n = (int)_init.size();
        info.assign((n << 2) + 1, Info());
        function<void(int, int, int)> build = [&](int p, int l, int r){
            if (l == r){
                info[p] = _init[l - 1];
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
 
    void modify(int p, int l, int r, int x, Info v){
        if (l == r){
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        if (x <= m){
            modify(2 * p, l, m, x, v);
        } 
        else{
            modify(2 * p + 1, m + 1, r, x, v);
        }
        pull(p);
    }
 
    void modify(int p, Info v){
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
        return query(2 * p, l, m, x, y) + query(2 * p + 1, m + 1, r, x, y);
    }
 
    Info query(int l, int r){
        return query(1, 1, n, l, r);
    }
};