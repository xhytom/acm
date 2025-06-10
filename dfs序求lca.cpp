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
// const ll MOD = 1e9+7;
int ksm(int x,int y){int ans=1;x%=MOD;while(y){if(y&1)ans=ans*x%MOD;x=x*x%MOD,y/=2;}return ans;}

const ll P1 = 999971, base1 = 101;
const ll P2 = 999973, base2 = 103;
const ll N = 200005;
//head


signed main()
{  
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    int n, m, s;
    cin >> n >> m >> s;
    vector<vector<int>> adj(n + 1);
    int u, v;
    for(int i = 1 ; i < n ; i++ )
    {
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    int idx = 0;
    vector<int> dfn(n + 5);
    vector st(__lg(n) + 2, vector<int> (n + 5));//****不能改成23****
    function<int(int,int)> get = [&](int x, int y) 
    {
        return dfn[x] < dfn[y] ? x : y;
    };
    function<void(int,int)> dfs = [&](int x, int fa) 
    {
        st[0][dfn[x] = ++idx] = fa;
        for(int y : adj[x]) if(y != fa) dfs(y, x); 
    };
    function<int(int,int)> lca = [&](int u, int v) 
    {
        if(u == v) return u;
        if((u = dfn[u]) > (v = dfn[v])) swap(u, v);
        int d = __lg(v - u++);
        return get(st[d][u], st[d][v - (1 << d) + 1]);
    };
    dfs(s, 0);
    for(int i = 1 ; i <= __lg(n) ; i++ )//****不能改成23****
    {
        for(int j = 1 ; j + (1 << i - 1) <= n  ; j++ ) // ****注意边界****
        {
            st[i][j] = get(st[i - 1][j], st[i - 1][j + (1 << i - 1)]);
        }
    }
    for(int i = 1 ; i <= m ; i++)
    {
        int u, v;
        cin >> u >> v;
        cout << lca(u, v) << "\n";
    }
    return 0;
}