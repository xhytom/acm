// sa[i]：排名为i的后缀的位置
// rk[i]：第i个位置开始的后缀的排名，作为基数排序的第一关键字 
// tp[i]：第二关键字中，排名为i的数的位置
// cnt[i]：有多少个元素排名为i
// s[i]：原输入数组 

// init: s[1..n], n = strlen(s + 1), m = SIGMA, makesa()
const int N = 1E6 + 5;
#define rep(i, s, t) for (int i = s; i <= t; ++i)
#define per(i, s, t) for (int i = t; i >= s; --i)

int m = 256;
char s[N];
int n, sa[N], rk[N], tp[N], cnt[N];
void init() {
    rep(i, 1, n) rk[i] = s[i], tp[i] = i;
}
void Qsort() {
    rep(i, 1, m) cnt[i] = 0;
    rep(i, 1, n) ++ cnt[rk[i]];
    rep(i, 1, m) cnt[i] += cnt[i - 1];
    per(i, 1, n) sa[cnt[rk[tp[i]]] --] = tp[i];
}
void get_sort() {
    for(int w = 1, p = 0; w <= n; m = p, p = 0, w <<= 1) {
        rep(i, n - w + 1, n) tp[++ p] = i;
        rep(i, 1, n) if(sa[i] > w) tp[++ p] = sa[i] - w;
        Qsort(), swap(rk, tp), p = rk[sa[1]] = 1;
        rep(i, 2, n) rk[sa[i]] = (tp[sa[i]] == tp[sa[i - 1]] && tp[sa[i] + w] == tp[sa[i - 1] + w]) ? p : ++ p;
        if(p == n) return;
    }
}
void makesa() {
    init();Qsort();get_sort();
}
int ht[N];
void makeht() {
    for(int k = 0, i = 1; i <= n; i++) {
        k = max(k - 1, 0);
        if(rk[i] == 1) continue;
        int j = sa[rk[i] - 1];
        while(s[i + k] == s[j + k]) k++;
        ht[rk[i]] = k;
    }
}


int st[25][N];
void makest() {
    rep(i, 1, n) st[0][i] = ht[i];
    rep(j, 1, 22) {
        for (int i = 1; i + (1 << j) - 1 <= n; i++) {
            st[j][i] = min(st[j - 1][i], st[j - 1][i + (1ll << j - 1)]);
        }
    }
}

int getLcp(int u, int v) {
    if(u == v) return n - u + 1;
    if(rk[u] > rk[v]) swap(u, v);
    int l = rk[u] + 1, r = rk[v];
    int len = __lg(r - l + 1);
    return min(st[len][l], st[len][r - (1 << len) + 1]);
}