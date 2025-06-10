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
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
void test() {cerr << "\n";}
template<typename T, typename... Args> 
void test(T x, Args... args) {cerr << x << " ";test(args...);}
const ll MOD = 998244353;
//const ll MOD = 1e9+7;
// ----------------------------✨✨✨---------------------------
const ll P1 = 999971, base1 = 101;
const ll P2 = 999973, base2 = 103;
const ll N = 30004;
//head

vector<int> e[N];
int fa[N], sz[N], dep[N], top[N], hs[N], dfn[N], idx, n, w[N];

struct Node{
    int mx = 0, sum = 0;
}tr[N << 2];

void build(int u, int L, int R)
{
    int mid = L + R >> 1;
    if(L == R)
    {
        tr[u].mx = 0;
        tr[u].sum = 0;
        return;
    }
    build(u << 1, L, mid);
    build(u << 1 | 1, mid + 1, R);
}

void pushup(int u, int L, int R)
{
    int mid = L + R >> 1;
    tr[u].sum = tr[u << 1].sum + tr[u << 1 | 1].sum;
    tr[u].mx = max(tr[u << 1].mx, tr[u << 1 | 1].mx);
}
void modify(int u, int L, int R, int x, int d)
{
    int mid = L + R >> 1;
    if(L == R)
    {
        tr[u].mx = d;
        tr[u].sum = d;
        return;
    }
    if(x <= mid)
    {
        modify(u << 1, L, mid, x, d);
    }else{
        modify(u << 1 | 1, mid + 1, R, x, d);
    }
    pushup(u, L, R);
}

int querySum(int u, int L, int R, int l, int r)
{
    int mid = L + R >> 1;
    int res = 0;
    if(l <= L && R <= r)
    {
        return tr[u].sum;
    }
    if(l <= mid)
    {
        res += querySum(u << 1, L, mid, l, r);
    }
    if(r > mid)
    {
        res += querySum(u << 1 | 1, mid + 1, R, l, r);
    }
    return res;
}

int queryMax(int u, int L, int R, int l, int r)
{
    int mid = L + R >> 1;
    int res = -1e18;
    if(l <= L && R <= r)
    {
        return tr[u].mx;
    }
    if(l <= mid)
    {
        res = max(res, queryMax(u << 1, L, mid, l, r));
    }
    if(r > mid)
    {
        res = max(res, queryMax(u << 1 | 1, mid + 1, R, l, r));
    }
    return res;
}

void dfs1(int x, int fx)
{
    sz[x] = 1;
    for(auto y : e[x])
    {
        if(y == fx) continue;
        dep[y] = dep[x] + 1;
        fa[y] = x;
        dfs1(y, x);
        sz[x] += sz[y];
        if(sz[y] > sz[hs[x]])
        {
            hs[x] = y;
        }
    }
}

void dfs2(int x, int root)
{
    dfn[x] = ++idx;
    top[x] = root;
    if(hs[x])
    {
        dfs2(hs[x], root);
    }
    for(auto y : e[x])
    {
        if(y == fa[x] || y == hs[x]) continue;
        dfs2(y, y);
    }
}

int lca(int x, int y)
{
    while(top[x] != top[y])
    {
        if(dep[top[x]] > dep[top[y]])
        {
            x = fa[top[x]];
        }else{
            y = fa[top[y]];
        }
    }
    return dep[x] > dep[y] ? y : x;
}

int query_sum(int x, int y)
{
    int ans = 0;
    while(top[x] != top[y])
    {
        if(dep[top[x]] > dep[top[y]])
        {
            ans += querySum(1, 1, n, dfn[top[x]], dfn[x]);
            x = fa[top[x]];
        }else{
            ans += querySum(1, 1, n, dfn[top[y]], dfn[y]);
            y = fa[top[y]];
        }
    }
    if(dep[x] < dep[y])
    {
        ans += querySum(1, 1, n, dfn[x], dfn[y]);
    }else{
        ans += querySum(1, 1, n, dfn[y], dfn[x]);
    }
    return ans;
}

int query_max(int x, int y)
{
    int ans = -1e18;
    while(top[x] != top[y])
    {
        if(dep[top[x]] > dep[top[y]])
        {
            ans = max(ans, queryMax(1, 1, n, dfn[top[x]], dfn[x]));
            x = fa[top[x]];
        }else{
            ans = max(ans,queryMax(1, 1, n, dfn[top[y]], dfn[y]));
            y = fa[top[y]];
        }
    }
    if(dep[x] < dep[y])
    {
        ans = max(ans,queryMax(1, 1, n, dfn[x], dfn[y]));
    }else{
        ans = max(ans,queryMax(1, 1, n, dfn[y], dfn[x]));
    }
    return ans;
}


signed main()
{  
#ifdef localfreopen
    freopen("1.in","r",stdin);
#endif
    fastio

    cin >> n;
    build(1, 1, n);
    for(int i = 1, u, v ; i < n ; i++ )
    {
        cin >> u >> v;
        e[u].pb(v);
        e[v].pb(u);
    }
    for(int i = 1 ; i <= n ; i++ )
    {
        cin >> w[i];
    }
    dfs1(1, 0);
    dfs2(1, 1);
    for(int i = 1 ; i <= n ; i++ )
    {
        modify(1, 1, n, dfn[i], w[i]);
    }
    int q;
    cin >> q;
    string s;
    int u, v, t;
    for(int i = 1 ; i <= q ; i++ )
    {
        cin >> s >> u >> v;
        if(s[0] == 'C')
        {
            modify(1, 1, n, dfn[u], v);
        }else if (s[1] == 'M')
        {
            cout << query_max(u, v) << "\n";
        }else{
            cout << query_sum(u, v) << "\n";
        }
    }
    return 0;
}