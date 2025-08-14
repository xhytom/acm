namespace INV {
    #define L(i, l, r) for (int i = l; i <= r; i++)
    #define R(i, r, l) for (int i = r; i >= l; i--)
    using ull = unsigned long long;
    using ll = long long;

    const int B = 1024, MOD = 998244353, N = (1 << 21) + 1;
    int K, pool[N * 2], *iv = pool + N;
    struct M {
        int ml, dis; 
    } mp[N];
    int inv(int w) {
        M d = mp[w >> 10];
        return (ull) iv[w * d.ml - d.dis] * (MOD + d.ml) % MOD;
    } 
    void init() {
        K = MOD / B;
        int Rlim = MOD - K;
        L(fb, 1, B) {
            int cur = 0, Add = fb * B;
            for(int p = 0; p <= K; ) {
                if(cur <= K) mp[p].ml = fb;
                else if(cur > Rlim) mp[p].ml = -fb;
                else {
                    int A = (Rlim - cur) / Add;
                    cur += A * Add, p += A;
                }
                cur += Add, ++p;
                if(cur >= MOD) cur -= MOD; 
            }
        }
        int count = 0;
        L(i, 1, K) 
            if(mp[i].ml > 0) mp[i].dis = (ll) mp[i].ml * i * B / MOD * MOD, ++count;
            else mp[i].dis = (ll) mp[i].ml * i * B / MOD * MOD - MOD, ++count;
        iv[1] = 1;
        L(i, 2, N - 1) iv[i] = (ll) iv[MOD % i] * (MOD - MOD / i) % MOD;
        L(i, 1, N - 1) iv[-i] = MOD - iv[i];
    }
}