#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int V = 1010;
const int E = 101000;
template<typename T>
struct MaxFlow
{
    int s, t, vtot;
    int head[V], etot;
    int dis[V], cur[V];
    struct edge
    {
        int v, nxt;
        T f;
    }e[E * 2];
    void addedge(int u, int v, T f)
    {
        e[etot] = {v, head[u], f}; head[u] = etot++;
        e[etot] = {u, head[v], 0}; head[v] = etot++;
    }
    bool bfs()
    {
        for(int i = 1 ; i <= vtot ; i++ )
        {
            dis[i] = 0;
            cur[i] = head[i];
        }
        queue<int> q;
        q.push(s); dis[s] = 1;
        while(!q.empty())
        {
            int u = q.front(); q.pop();
            for(int i = head[u] ; ~i ; i = e[i].nxt)
            {
                if(e[i].f && !dis[e[i].v])
                {
                    int v = e[i].v;
                    dis[v] = dis[u] + 1;
                    if(v == t) return true;
                    q.push(v);
                }
            }
        }
        return false;
    }
    T dfs(int u, T m)
    {
        if(u == t) return m;
        T flow = 0;
        for(int i = cur[u]; ~i ; cur[u] = i = e[i].nxt)
        {
            if(e[i].f && dis[e[i].v] == dis[u] + 1)
            {
                T f = dfs(e[i].v, min(m, e[i].f));
                e[i].f -= f;
                e[i ^ 1].f += f;
                m -= f;
                flow += f;
                if(!m) break;
            }
        }
        if(!flow) dis[u] = -1;
        return flow;
    }
    T dinic()
    {
        T flow = 0;
        while(bfs()) flow += dfs(s, numeric_limits<T>::max());
        return flow;
    }
    void init(int s_, int t_, int vtot_ )
    {
        s = s_;
        t = t_;
        vtot = vtot_;
        etot = 0;
        for(int i = 1 ; i <= vtot ; i++ )
        {
            head[i] = -1;
        }
    } 
};

MaxFlow<ll> g;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    int n, m, s, t;
    cin >> n >> m >> s >> t;
    g.init(s, t, n);
    for(ll i = 1, u, v, w ; i <= m ; i++ )
    {
        cin >> u >> v >> w;
        g.addedge(u, v, w);
    }
    cout << g.dinic() << "\n";
}
