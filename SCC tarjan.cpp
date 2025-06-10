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
const ll N = 100005;
//head

vector<int> e[N];
int dfn[N], ins[N], low[N], bel[N], idx, cnt;
stack<int> st;
vector<vector<int> > scc;


void dfs(int u)
{
    dfn[u] = low[u] = ++idx;
    ins[u] = true;
    st.push(u);
    for(auto v : e[u])
    {
        if(!dfn[v])
        {
            dfs(v);
            low[u] = min(low[u], low[v]);
        }else{
            if(ins[v]) low[u] = min(low[u], dfn[v]);
        }
    }
    if(dfn[u] == low[u])
    {
        vector<int> c;
        ++cnt;
        while(true)
        {
            int v = st.top();
            c.push_back(v);
            ins[v] = false;
            bel[v] = cnt;
            st.pop();
            //cout << v << " ";
            if(v == u) break;
        }
        //cout << endl;
        sort(c.begin(), c.end());
        scc.push_back(c);
    }

}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    int n, m;
    cin >> n >> m;
    for(int i = 1 ; i <= m ; i++ )
    {
        int x, y;
        cin >> x >> y;
        e[x].push_back(y);
    }
    for(int i = 1 ; i <= n ; i++ )
    {
        if(!dfn[i])
        {
            dfs(i);
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