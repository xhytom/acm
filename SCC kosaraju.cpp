#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
#define debug(x) cerr << #x << " = " << (x) << endl;
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
void test() {cerr << "\n";}
template<typename T, typename... Args> 
void test(T x, Args... args) {cerr << x << " ";test(args...);}
const ll MOD = 998244353;
//const ll MOD = 1e9+7;
const ll N = 200005;
//head

int vis[N], n, m;
vector<int> out, c, e[N], erev[N];
int sz[N];
int bel[N], cnt;
vector<vector<int> >scc;

void dfs1(int u)
{
    vis[u] = 1;
    for(auto v : e[u])
    {
        if(!vis[v]) dfs1(v);
    }
    out.push_back(u);
}

void dfs2(int u, int cnt)
{

    vis[u] = 1;
    for(auto v : erev[u])
    {
        if(!vis[v]) dfs2(v, cnt);
    }
    bel[u] = cnt;
    sz[cnt]++;
    c.push_back(u);
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    int n, m, x, y;
    cin >> n >> m;
    for(int i = 1 ; i <= m ; i++ )
    {
        cin >> x >> y;
        e[x].push_back(y);
        erev[y].push_back(x);
    }
    memset(vis, 0, sizeof(vis));
    for(int i = 1 ; i <= n ; i++ )
    {
        if(!vis[i])
        {
            dfs1(i);
        }
    }
    reverse(out.begin(), out.end());
    memset(vis, 0, sizeof(vis));
    for(auto u : out)
    {
        if(!vis[u])
        {
            c.clear();
            dfs2(u, ++cnt);
            sort(c.begin(), c.end());
            scc.push_back(c);
        }
        
    }
    sort(scc.begin(), scc.end());
    for(auto c : scc)
    {
        for(auto x : c)
        {
            cout << x << " ";
        }
        cout << "\n";
    }
    return 0;
}