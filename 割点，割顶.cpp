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
const ll N = 20005;
//head

int n, m;
int dfn[N], idx, low[N];
bool vis[N], cut[N];
vector<int> e[N];
int cnt;

void dfs(int u, int root)
{
    vis[u] = 1;
    dfn[u] = ++idx;
    low[u] = idx;
    int child = 0;
    for(auto v : e[u])
    {
        if(!vis[v])
        {
            dfs(v, root);
            low[u] = min(low[u], low[v]);
            if(low[v] >= dfn[u] && u != root)
            {
                cut[u] = 1;
            }
            if(u == root)
            {
                child++;
            }
        }
        low[u] = min(low[u], dfn[v]);
    }
    if(child >= 2 && u == root)
    {
        cut[u] = 1;
    }
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    cin >> n >> m;
    rep(i, 1, m + 1)
    {
        int x, y;
        cin >> x >> y;
        e[x].push_back(y);
        e[y].push_back(x);
    }
    rep(i, 1, n + 1)
    {
        if(!vis[i])
        {
            dfs(i, i);
        }
    }
    cout << accumulate(cut + 1, cut + n + 1, 0ll) << "\n";
    rep(i, 1, n + 1)
    {
        if(cut[i])
        {
            cout << i << " ";
        }
    }
    return 0;
}