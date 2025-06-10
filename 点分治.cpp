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
    int n, k, ans = 0;
    cin >> n >> k;
    ans = n + 1;
    vector<vector<pair<int,int>>> adj(n + 1);
    vector<int> sz(n + 1, 0), maxsz(n + 1, 0), del(n + 1, 0);
    vector<int> mark(k + 1, 0), c(k + 1, 0);
    int T = 1;
    int u, v, w;
    for(int i = 1 ; i < n ; i++ )
    {
        cin >> u >> v >> w;
        u++;
        v++;
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }
    function<void(int, int)> solve = [&](int x, int s)
    {
        T++;
        int mxs = s + 1, root = -1;
        function<void(int, int)> dfs1 = [&](int x, int fx)
        {
            sz[x] = 1;
            maxsz[x] = 0;
            for(auto [y, w] : adj[x])
            {
                if(del[y] || y == fx) continue;
                dfs1(y, x);
                sz[x] += sz[y];
                maxsz[x] = max(maxsz[x], sz[y]);
            }
            maxsz[x] = max(maxsz[x], s - sz[x]);
            if(maxsz[x] < mxs)
            {
                mxs = maxsz[x], root = x;
            }
        };
        dfs1(x, -1);
        // cout << root << endl;
        /////////////////////////////////
        mark[0] = T;
        c[0] = 0;
        for(auto [y, w] : adj[root])
        {
            if(del[y]) continue;
            vector<pair<int, int>> self;
            function<void(int, int, int, int)> dfs2 = [&](int x, int fx, int dis, int dep)
            {
                self.emplace_back(dis, dep);
                for(auto [y, w] : adj[x])
                {
                    if(del[y] || y == fx) continue;
                    dfs2(y, x, dis + w, dep + 1);
                }
            };
            dfs2(y, root, w, 1);
            for(auto [dis, dep] : self)
            {
                if(k - dis >= 0 && mark[k - dis] == T)
                {
                    ans = min(ans, c[k - dis] + dep);
                }
            }
            for(auto [dis, dep] : self)
            {
                if(dis > k) continue;
                if(mark[dis] == T)
                {
                    c[dis] = min(c[dis], dep);
                }else{
                    c[dis] = dep;
                    mark[dis] = T;
                }
            }
        }
        /////////////////////////////////
        del[root] = 1;
        for(auto [y, w] : adj[root])
        {
            if(del[y]) continue;
            solve(y, sz[y]);
        }
    };
    solve(1, n);
    cout << (ans > n ? -1 : ans) << "\n";
    return 0;
}