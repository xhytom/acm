/*
 
_/      _/    _/      _/    _/      _/   _/_/_/_/_/     _/_/       _/      _/ 
 _/    _/     _/      _/     _/    _/        _/       _/    _/     _/      _/            
  _/  _/      _/      _/      _/  _/         _/      _/      _/    _/_/  _/_/         
   _/_/       _/_/_/_/_/        _/           _/      _/      _/    _/  _/  _/          
  _/  _/      _/      _/        _/           _/      _/      _/    _/      _/          
 _/    _/     _/      _/        _/           _/       _/    _/     _/      _/          
_/      _/    _/      _/        _/           _/         _/_/       _/      _/       
 
*/
#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
#define debug(x) cerr << #x << " = " << (x) << endl;
#define int long long
#define pb push_back
#define eb emplace_back
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937_64 mrand(chrono::steady_clock().now().time_since_epoch().count());
int rnd(int l,int r){ return mrand() % (r - l + 1) + l;}
void test() {cerr << "\n";}
template<typename T, typename... Args> 
void test(T x, Args... args) {cerr << x << " ";test(args...);}
const ll MOD = 998244353;
// const ll MOD = 1e9+7;
ll ksm(ll x,ll y){ll ans=1;x%=MOD;while(y){if(y&1)ans=ans*x%MOD;x=x*x%MOD,y/=2;}return ans;}

const ll P1 = 999971, base1 = 101;
const ll P2 = 999973, base2 = 103;
const ll N = 200005;
//head

struct DSU {
    std::vector<int> f, siz;
    DSU(int n) : f(n), siz(n, 1) { std::iota(f.begin(), f.end(), 0); }
    int leader(int x) {
        while (x != f[x]) x = f[x] = f[f[x]];
        return x;
    }
    bool same(int x, int y) { return leader(x) == leader(y); }
    bool merge(int x, int y) {
        x = leader(x);
        y = leader(y);
        if (x == y) return false;
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
    int size(int x) { return siz[leader(x)]; }
};

struct CutEdges {
    int n;
    int idx = 0;
    int iddx = 1;
    vector<int> low, dfn, fa, head, to, nxt;
    vector<pair<int,int>> bridge;
    CutEdges(int n, int m) : low(n + 1), dfn(n + 1), fa(n + 1), adj(n + 1) 
    head(2 * m + 4), nxt(2 * m + 4), to(2 * m + 4) {
        this->n = n;
    }
    void addEdge(int x, int y) {
        nxt[++iddx] = head[x];
        head[x] = iddx;
        to[iddx] = y;
    }
    vector<pair<int, int>> work() {
        for (int i = 1; i <= n; i++) {
            if (!dfn[i]) tarjan(i, 0);
        }
        return bridge;
    }
    void tarjan(int x, int e_in)
{
    dfn[x] = low[x] = ++idx;
    for(int i = head[x]; i; i = nxt[i]) {
        int y = to[i];
        if(!dfn[y]) {
            tarjan(y, i);
            if(dfn[x] < low[y]) {
                b[i] = b[i ^ 1] = 1;
            }
            low[x] = min(low[x], low[y]);
        } else if (i != (e_in ^ 1)) {
            low[x] = min(low[x], dfn[y]);
        }
  
};

signed main()
{  
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    int n, m;
    std::cin >> n >> m;
    CutEdges g(n, m);
    vector<pair<int,int>> edges(m);
    for (int i = 0; i < m; i++) {
        int u, v;
        std::cin >> u >> v;
        edges[i] = {u, v};
        g.addEdge(u, v);
        g.addEdge(v, u);
    }
    auto bridge = g.work();
    debug(bridge.size());
    map<pair<int,int>, int> is;
    for (auto [x, y] : bridge) {
        is[{x, y}] = 1;
        is[{y, x}] = 1;
    }
    DSU dsu(n + 1);
    for (auto [x, y] : edges) {
        if (is[{x, y}] == 0) {
            dsu.merge(x, y);
        }
    }
    vector<vector<int>> bcc(n + 1);
    int ans = 0;
    for (int i = 1; i <= n; i++) {
        bcc[dsu.leader(i)].push_back(i);
        if (dsu.leader(i) == i) ans++;
    }
    std:;cout << ans << "\n";
    for (int i = 1; i <= n; i++) {
        if (bcc[i].size() >= 1) {
            cout << bcc[i].size() << " ";
            for (auto x : bcc[i]) {
                cout << x << " ";
            }
            cout << "\n";
        }
    }
    return 0;
}