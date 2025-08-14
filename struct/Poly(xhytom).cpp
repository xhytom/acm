const i32 MOD = 998244353;
using Z = MInt<MOD>;

std::vector<Z> roots{0, 1};
std::vector<int> rev;

void dft(std::vector<Z> &a){
    i32 n = a.size();
    if(rev.size() != n){
        rev.resize(n);
        i32 k = __builtin_ctz(n) - 1;
        for (i32 i = 0; i < n; i++) rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
    for (i32 i = 0; i < n; i++) if (i < rev[i]) std::swap(a[i], a[rev[i]]);
    if (roots.size() < n) {
        i32 k = __builtin_ctz(roots.size());
        roots.resize(n);
        while((1 << k) < n){
            Z e = Z(3) ^ ((MOD - 1) >> (k + 1));
            for(i32 i = 1 << (k - 1); i < 1 << k; i++) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = roots[i] * e;
            }
            k++;
        }
    }
    for (i32 k = 1; k < n; k *= 2) {
        for (i32 i = 0; i < n; i += 2 * k) {
            for (i32 j = 0; j < k; j++) {
                Z u = a[i + j], v = a[i + j + k] * roots[k + j];
                a[i + j] = u + v;
                a[i + j + k] = u - v;
            }
        }
    }
}
void idft(std::vector<Z> &a) {
    std::reverse(a.begin() + 1, a.end());
    dft(a);
    i32 n = a.size();
    Z inv = Z((1 - MOD) / n + MOD);
    for (i32 i = 0; i < n; i++) a[i] = a[i] * inv;
}
struct Poly : std::vector<Z> {

    Poly() {}
    Poly(const int &n) : std::vector<Z>(n) {}
    Poly(const std::vector<Z> &a) : std::vector<Z>(a) {}
    template<class It>
    explicit Poly(It first, It last) : std::vector<Z>(first, last) {}

    friend Poly operator*(Poly a, Poly b){
        if (a.empty() || b.empty())return Poly();
        if (a.size() > b.size()) std::swap(a, b);
        if (a.size() < 4) {
            Poly c(a.size() + b.size() - 1);
            for(i32 i = 0; i < a.size(); i++) {
                for(i32 j = 0; j < b.size(); j++) {
                    c[i + j] = c[i + j] + a[i] * b[j];
                }
            }
            return c;
        }
        i32 sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot) sz *= 2;
        a.resize(sz); b.resize(sz);
        dft(a); dft(b);
        for(i32 i = 0; i < sz; i++) a[i] = a[i] * b[i];
        idft(a);
        a.resize(tot);
        return a;
    }
};