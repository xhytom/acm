ll prim[50000005], sum[50000005], d[50000005], len;
bool vis[50000005];

inline void sieve(int x) {
    for(int i = 2;i <= x;i ++) {
        if(! vis[i]) {
            prim[++ len] = i;
            d[i] = 2;
            sum[i] = 1;
        }
        for(int j = 1;j <= len && i * prim[j] <= x;j ++) {
            vis[i * prim[j]] = 1;
            if(i % prim[j] == 0) {
                sum[i * prim[j]] = sum[i] + 1;
                d[i * prim[j]] = d[i] / (sum[i] + 1) * (sum[i] + 2);
                break;
            }
            sum[i * prim[j]] = 1;
            d[i * prim[j]] = d[i] * 2;
        }
    }
}
