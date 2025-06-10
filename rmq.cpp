template<typename T>
struct Rmq {
    std::vector<std::vector<T>> f;
    Rmq (std::vector<T> &q) {
        int n = q.size();
        f.assign(std::__lg(n) + 1, std::vector<T>(n));
        for (int i = 0; i < n; i++) {
            f[0][i] = q[i];
        }
        for (int i = 1; i <= std::__lg(n); i++) {
            for(int j = 0 ; j + (1 << i) - 1 < n ; j++ ) {
                f[i][j] = f[i - 1][j] > f[i - 1][j + (1 << i - 1)] ? f[i - 1][j] : f[i - 1][j + (1 << i - 1)];
            }
        }
    }

    T query(int l, int r) {
        int len = std::__lg(r - l + 1);
        return f[len][l] > f[len][r - (1 << len) + 1] ? f[len][l] : f[len][r - (1 << len) + 1];
    }

};