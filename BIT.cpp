template<class T>
struct BIT {
    int size = 1;
    std::vector<T> c;
    BIT (int x) {
        size = x + 5;
        c.resize(x + 5);
    }
    void init(int x) {
        size = x + 5;
        c.resize(x + 5);
    }
    void clear() {
        for(int i = 0; i < size; i++) {
            c[i] = 0;
        }
    }
    void change(int x, T y) {
        for(; x < size ; x += x & (-x)) {
            c[x] += y;
        }
    }
    T query(int x) {
        T s = 0;
        for(;x ;x -= x & (-x)) {
            s += c[x];
        }
        return s;
    }
    T query(int l, int r) {
        if (l == 0) return query(r);
        return query(r) - query(l - 1);
    }
    int kth(int k) {
        int sum = 0, x = 0;
        for (int i = log2(size); ~i; --i) {
            x += 1 << i;
            if (x >= size || sum + c[x] >= k)
                x -= 1 << i;
            else
                sum += c[x];
        }
        return x + 1;
    }
};