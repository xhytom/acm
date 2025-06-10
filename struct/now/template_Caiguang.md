[TOC] 

# 一、杂项

## compare

```cpp
:loop
 data.exe > 1.in
 my.exe <1.in >my.out
 std.exe <1.in >std.out
 fc my.out std.out
if not errorlevel 1 goto loop
pause
goto loop
```

## 怪分组背包

```cpp
void solve() {
	int n, v;
	cin >> n >> v;
	vector<int> s(n + 1);
	vector V(n + 1, vector<int>()), w(n + 1, vector<int>());
	for(int i = 1; i <= n; i++) {
		cin >> s[i];
		int x;
		for(int j = 1; j <= s[i]; j++) {
			cin >> x;
			V[i].push_back(x);
			cin >> x;
			w[i].push_back(x);
		}
	}
	vector<int> dp(v + 1);
	int now = 0;
	for(int i = 1; i <= n; i++) {
		now ^= 1;
		for(int k = v; k >= 0; k--) {
			for(int j = 0; j < s[i]; j++) {
				if(k >= V[i][j])
				dp[k] = max(dp[k], dp[k - V[i][j]] + w[i][j]);
			}
		}
	}
	int ma = 0;
	cout << dp[v] << '\n';
}

```



## __int128

```c++
__int128 read()
{
    __int128 f=1,w=0;
    char ch=getchar();
    while(ch<'0'||ch>'9')
    {
        if(ch=='-')
        f=-1;
        ch=getchar();
    }
    while(ch<='9'&&ch>='0')
    {
        w=w*10+ch-'0';
        ch=getchar();
    }
    return f*w;
}

void print(__int128 x)
{
    if(x<0)
    {
        putchar('-');
        x=-x;
    }
    if(x>9)print(x/10);
    putchar(x%10+'0');
}
```

## exgcd

```cpp
using i64 = long long;
 
i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    i64 g = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return g;
}
```

## DSU

```cpp
struct DSU {
    std::vector<int> f, siz;
     
    DSU() {}
    DSU(int n) {
        init(n);
    }
     
    void init(int n) {
        f.resize(n);
        std::iota(f.begin(), f.end(), 0);
        siz.assign(n, 1);
    }
     
    int find(int x) {
        while (x != f[x]) {
            x = f[x] = f[f[x]];
        }
        return x;
    }
     
    bool same(int x, int y) {
        return find(x) == find(y);
    }
     
    bool merge(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            return false;
        }
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
     
    int size(int x) {
        return siz[find(x)];
    }
};
```

## simpson

```cpp
double f(double x) {
    return ;
}
double simpson(double l, double r) {
    return (f(l) + 4 * f((l + r) / 2) + f(r)) * (r - l) / 6;
}
double integral(double l, double r, double eps, double st) {
    double mid = (l + r) / 2;
    double sl = simpson(l, mid);
    double sr = simpson(mid, r);
    if (std::abs(sl + sr - st) <= 15 * eps)
        return sl + sr + (sl + sr - st) / 15;
    return integral(l, mid, eps / 2, sl) + integral(mid, r, eps / 2, sr);
}
double integral(double l, double r) {
    return integral(l, r, EPS, simpson(l, r));
}
```

## linear_basis

```cpp
struct linearBasis {
    vector<int> a;
    int sz;
    linearBasis() : sz(0), a(70) {}
    void insert(int x) {
        for(int i = 60; i >= 0; i--) {
            if(!(x&(1ll<<i))) continue;
            if(a[i]) {
                x ^= a[i];
            } else {
                a[i] = x;
                sz++;
                return;
            }
        }
    }
};
```

## poly

```cpp
#include <bits/stdc++.h>
using i64 = long long;
using u64 = unsigned long long;
using u32 = unsigned;
constexpr int P = 998244353;
std::vector<int> rev, roots{0, 1};
int power(int a, int b) {
    int res = 1;
    for (; b; b >>= 1, a = 1ll * a * a % P)
        if (b & 1)
            res = 1ll * res * a % P;
    return res;
}
void dft(std::vector<int> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; ++i)
            rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
    for (int i = 0; i < n; ++i)
        if (rev[i] < i)
            std::swap(a[i], a[rev[i]]);
    if (int(roots.size()) < n) {
        int k = __builtin_ctz(roots.size());
        roots.resize(n);
        while ((1 << k) < n) {
            int e = power(3, (P - 1) >> (k + 1));
            for (int i = 1 << (k - 1); i < (1 << k); ++i) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = 1ll * roots[i] * e % P;
            }
            ++k;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; ++j) {
                int u = a[i + j];
                int v = 1ll * a[i + j + k] * roots[k + j] % P;
                int x = u + v;
                if (x >= P)
                    x -= P;
                a[i + j] = x;
                x = u - v;
                if (x < 0)
                    x += P;
                a[i + j + k] = x;
            }
        }
    }
}
void idft(std::vector<int> &a) {
    int n = a.size();
    std::reverse(a.begin() + 1, a.end());
    dft(a);
    int inv = power(n, P - 2);
    for (int i = 0; i < n; ++i)
        a[i] = 1ll * a[i] * inv % P;
}
struct Poly {
    std::vector<int> a;
    Poly() {}
    Poly(int a0) {
        if (a0)
            a = {a0};
    }
    Poly(const std::vector<int> &a1) : a(a1) {
        while (!a.empty() && !a.back())
            a.pop_back();
    }
    int size() const {
        return a.size();
    }
    int operator[](int idx) const {
        if (idx < 0 || idx >= size())
            return 0;
        return a[idx];
    }
    Poly mulxk(int k) const {
        auto b = a;
        b.insert(b.begin(), k, 0);
        return Poly(b);
    }
    Poly modxk(int k) const {
        k = std::min(k, size());
        return Poly(std::vector<int>(a.begin(), a.begin() + k));
    }
    Poly divxk(int k) const {
        if (size() <= k)
            return Poly();
        return Poly(std::vector<int>(a.begin() + k, a.end()));
    }
    friend Poly operator+(const Poly a, const Poly &b) {
        std::vector<int> res(std::max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); ++i) {
            res[i] = a[i] + b[i];
            if (res[i] >= P)
                res[i] -= P;
        }
        return Poly(res);
    }
    friend Poly operator-(const Poly a, const Poly &b) {
        std::vector<int> res(std::max(a.size(), b.size()));
        for (int i = 0; i < int(res.size()); ++i) {
            res[i] = a[i] - b[i];
            if (res[i] < 0)
                res[i] += P;
        }
        return Poly(res);
    }
    friend Poly operator*(Poly a, Poly b) {
        int sz = 1, tot = a.size() + b.size() - 1;
        while (sz < tot)
            sz *= 2;
        a.a.resize(sz);
        b.a.resize(sz);
        dft(a.a);
        dft(b.a);
        for (int i = 0; i < sz; ++i)
            a.a[i] = 1ll * a[i] * b[i] % P;
        idft(a.a);
        return Poly(a.a);
    }
    Poly &operator+=(Poly b) {
        return (*this) = (*this) + b;
    }
    Poly &operator-=(Poly b) {
        return (*this) = (*this) - b;
    }
    Poly &operator*=(Poly b) {
        return (*this) = (*this) * b;
    }
    Poly deriv() const {
        if (a.empty())
            return Poly();
        std::vector<int> res(size() - 1);
        for (int i = 0; i < size() - 1; ++i)
            res[i] = 1ll * (i + 1) * a[i + 1] % P;
        return Poly(res);
    }
    Poly integr() const {
        if (a.empty())
            return Poly();
        std::vector<int> res(size() + 1);
        for (int i = 0; i < size(); ++i)
            res[i + 1] = 1ll * a[i] * power(i + 1, P - 2) % P;
        return Poly(res);
    }
    Poly inv(int m) const {
        Poly x(power(a[0], P - 2));
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (2 - modxk(k) * x)).modxk(k);
        }
        return x.modxk(m);
    }
    Poly log(int m) const {
        return (deriv() * inv(m)).integr().modxk(m);
    }
    Poly exp(int m) const {
        Poly x(1);
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x * (1 - x.log(k) + modxk(k))).modxk(k);
        }
        return x.modxk(m);
    }
    Poly sqrt(int m) const {
        Poly x(1);
        int k = 1;
        while (k < m) {
            k *= 2;
            x = (x + (modxk(k) * x.inv(k)).modxk(k)) * ((P + 1) / 2);
        }
        return x.modxk(m);
    }
    Poly mulT(Poly b) const {
        if (b.size() == 0)
            return Poly();
        int n = b.size();
        std::reverse(b.a.begin(), b.a.end());
        return ((*this) * b).divxk(n - 1);
    }
    std::vector<int> eval(std::vector<int> x) const {
        if (size() == 0)
            return std::vector<int>(x.size(), 0);
        const int n = std::max(int(x.size()), size());
        std::vector<Poly> q(4 * n);
        std::vector<int> ans(x.size());
        x.resize(n);
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                q[p] = std::vector<int>{1, (P - x[l]) % P};
            } else {
                int m = (l + r) / 2;
                build(2 * p, l, m);
                build(2 * p + 1, m, r);
                q[p] = q[2 * p] * q[2 * p + 1];
            }
        };
        build(1, 0, n);
        std::function<void(int, int, int, const Poly &)> work = [&](int p, int l, int r, const Poly &num) {
            if (r - l == 1) {
                if (l < int(ans.size()))
                    ans[l] = num[0];
            } else {
                int m = (l + r) / 2;
                work(2 * p, l, m, num.mulT(q[2 * p + 1]).modxk(m - l));
                work(2 * p + 1, m, r, num.mulT(q[2 * p]).modxk(r - m));
            }
        };
        work(1, 0, n, mulT(q[1].inv(n)));
        return ans;
    }
};

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    return 0;
}
```

# 二、网络流

## flow

```cpp
#include <bits/stdc++.h>

using i64 = long long;

template<typename T>
class FlowGraph {
public:
    const int V = 10100;
    const int E = 101000;
    using ll = long long;
    int s, t, vtot;
    int head[V], etot;
    int dis[V], cur[V];
    struct edge {
        int v, nxt;
        T f;
    } e[E];
    void addedge(int u, int v, T f) {
        e[etot] = {v, head[u], f}; head[u] = etot++;
        e[etot] = {u, head[v], 0}; head[v] = etot++;
    }
    bool bfs() {
        for(int i = 1 ; i <= vtot ; i++ ) {
            dis[i] = 0;
            cur[i] = head[i];
        }
        std::queue<int> q;
        q.push(s); dis[s] = 1;
        while(!q.empty()) {
            int u = q.front(); q.pop();
            for(int i = head[u] ; ~i ; i = e[i].nxt) {
                if(e[i].f && !dis[e[i].v]) {
                    int v = e[i].v;
                    dis[v] = dis[u] + 1;
                    if(v == t) return true;
                    q.push(v);
                }
            }
        }
        return false;
    }
    T dfs(int u, T m) {
        if(u == t) return m;
        T flow = 0;
        for(int i = cur[u]; ~i ; cur[u] = i = e[i].nxt) {
            if(e[i].f && dis[e[i].v] == dis[u] + 1) {
                T f = dfs(e[i].v, std::min(m, e[i].f));
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
    T dinic() {
        T flow = 0;
        while(bfs()) flow += dfs(s, std::numeric_limits<T>::max());
        return flow;
    }
    void init(int s_, int t_, int vtot_ ) {
        s = s_;
        t = t_;
        vtot = vtot_;
        etot = 0;
        for(int i = 1 ; i <= vtot + 100 ; i++ ) {
            head[i] = -1;
        }
    } 
};

//***记得每次init,

FlowGraph<int> g;

```

## costflow

```cpp
#include <bits/stdc++.h>

using namespace std;

const int V = 2010;
const int E = 20100;
// #define int double
using ll = long long;

template<typename T>
struct MaxFlow {
    int s, t, vtot;
    int head[V], etot, cur[V];
    int pre[V];
    bool vis[V];
    T dis[V], cost, flow;

    struct edge {
        int v, nxt;
        T f, c;
    }e[E * 2];

    void addedge(int u, int v, T f, T c, T f2 = 0)
    {
        e[etot] = {v, head[u], f, c}; head[u] = etot++;
        e[etot] = {u, head[v], f2, -c}; head[v] = etot++;
    }

    bool spfa() {
        T inf = numeric_limits<T>::max() / 2;
        for(int i = 1; i <= vtot; i++) {
            dis[i] = inf;
            vis[i] = false;
            pre[i] = -1;
            cur[i] = head[i];
        }
        dis[s] = 0;
        vis[s] = true;
        queue<int> q;
        q.push(s);
        while(!q.empty()) {
            int u = q.front();
            for(int i = head[u]; ~i; i = e[i].nxt) {
                int v = e[i].v;
                if(e[i].f && dis[v] > dis[u] + e[i].c) {
                    dis[v] = dis[u] + e[i].c;
                    pre[v] = i;
                    if(!vis[v]) {
                        vis[v] = 1;
                        q.push(v);
                    } 
                }
            }
            q.pop();
            vis[u] = false;
        }
        return dis[t] < inf;
    }

    void augment() {
        int u = t;
        T f = numeric_limits<T>::max();
        while(~pre[u]) {
            f = min(f, e[pre[u]].f);
            u = e[pre[u] ^ 1].v;
        }
        flow += f;
        cost += f * dis[t];
        u = t;
        while(~pre[u]) {
            e[pre[u]].f -= f;
            e[pre[u] ^ 1].f += f;
            u = e[pre[u] ^ 1].v;
        }
    }

    pair<T, T> sol() {
        flow = cost = 0;
        while(spfa()) {
            augment();
        }
        return {flow, cost};
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

//***记得每次init,

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.tie(nullptr);
    int n, m;
    cin >> n >> m;
    int num = 0;
    vector<vector<int>> a(m + 1, vector<int> (n + m + 1)), cnt(m + 1, vector<int> (n + m + 1));
    for(int i = 1, j = 0; i <= m; i++, j++) {
        for(int k = 1; k <= n + j; k++) {
            cin >> a[i][k];
            a[i][k] *= -1;
            cnt[i][k] = ++num;
        }
    }
    int s = num * 2 + 1, t = s + 1;
    MaxFlow<double> g1, g2, g3;
    g1.init(s, t, t);
    g2.init(s, t, t);
    g3.init(s, t, t);
    for(int i = 1, j = 0; i <= m; i++, j++) {
        for(int k = 1; k <= n + j; k++) {
            g1.addedge(cnt[i][k], cnt[i][k] + num, 1, a[i][k]);
            g2.addedge(cnt[i][k], cnt[i][k] + num, n, a[i][k]);
            g3.addedge(cnt[i][k], cnt[i][k] + num, n, a[i][k]);
            if(i == 1) { 
                g1.addedge(s, cnt[i][k], 1, 0);
                g2.addedge(s, cnt[i][k], 1, 0);
                g3.addedge(s, cnt[i][k], 1, 0);
            } 
            if(i == m) {
                g1.addedge(cnt[i][k] + num, t, 1, 0);
                g2.addedge(cnt[i][k] + num, t, n, 0);
                g3.addedge(cnt[i][k] + num, t, n, 0);
            }
                if(k > 1 && i > 1) {
                    g1.addedge(cnt[i - 1][k - 1] + num, cnt[i][k], 1, 0);
                    g2.addedge(cnt[i - 1][k - 1] + num, cnt[i][k], 1, 0);
                    g3.addedge(cnt[i - 1][k - 1] + num, cnt[i][k], n, 0);
                }
                if(k != n + j && i > 1) {
                    g1.addedge(cnt[i - 1][k] + num, cnt[i][k], 1, 0);
                    g2.addedge(cnt[i - 1][k] + num, cnt[i][k], 1, 0);
                    g3.addedge(cnt[i - 1][k] + num, cnt[i][k], n, 0);
                }
        }
    }
    cout << fabs(g1.sol().second) << '\n' << fabs(g2.sol().second) << '\n' << fabs(g3.sol().second);
}
```



## jly_flow

```cpp
using i64 = long long;
template<class T>
struct MaxFlow {
    struct _Edge {
        int to;
        T cap;
        _Edge(int to, T cap) : to(to), cap(cap) {}
    };
    
    int n;
    std::vector<_Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<int> cur, h;
    
    MaxFlow() {}
    MaxFlow(int n) {
        init(n);
    }
    
    void init(int n) {
        this->n = n;
        e.clear();
        g.assign(n, {});
        cur.resize(n);
        h.resize(n);
    }
    
    bool bfs(int s, int t) {
        h.assign(n, -1);
        std::queue<int> que;
        h[s] = 0;
        que.push(s);
        while (!que.empty()) {
            const int u = que.front();
            que.pop();
            for (int i : g[u]) {
                auto [v, c] = e[i];
                if (c > 0 && h[v] == -1) {
                    h[v] = h[u] + 1;
                    if (v == t) {
                        return true;
                    }
                    que.push(v);
                }
            }
        }
        return false;
    }
    
    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        auto r = f;
        for (int &i = cur[u]; i < int(g[u].size()); ++i) {
            const int j = g[u][i];
            auto [v, c] = e[j];
            if (c > 0 && h[v] == h[u] + 1) {
                auto a = dfs(v, t, std::min(r, c));
                e[j].cap -= a;
                e[j ^ 1].cap += a;
                r -= a;
                if (r == 0) {
                    return f;
                }
            }
        }
        return f - r;
    }
    void addEdge(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }
    T flow(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, std::numeric_limits<T>::max());
        }
        return ans;
    }
    
    std::vector<bool> minCut() {
        std::vector<bool> c(n);
        for (int i = 0; i < n; i++) {
            c[i] = (h[i] != -1);
        }
        return c;
    }
    
    struct Edge {
        int from;
        int to;
        T cap;
        T flow;
    };
    std::vector<Edge> edges() {
        std::vector<Edge> a;
        for (int i = 0; i < e.size(); i += 2) {
            Edge x;
            x.from = e[i + 1].to;
            x.to = e[i].to;
            x.cap = e[i].cap + e[i + 1].cap;
            x.flow = e[i + 1].cap;
            a.push_back(x);
        }
        return a;
    }
};
```

```cpp
#include <bits/stdc++.h>

using i64 = long long;

template <class T>
struct Maxflow {
	struct _Edge {
		int to;
		T cap;
		_Edge(int to, T cap) : to(to), cap(cap) {}
	};

	int n;
	std::vector<_Edge> e;
	std::vector<std::vector<int>> g;
	std::vector<int> cur, h;

	Maxflow() {}

	Maxflow(int n) {
		init(n);
	}

	void init(int n) {
		this -> n = n;
		e.clear();
		g.assign(n + 1, {});
		cur.resize(n + 1);
		h.resize(n + 1);
	}

	bool bfs(int s, int t) {
		h.assign(n + 1, -1);
		std::queue<int> que;
		h[s] = 0;
		que.push(s);
		while(que.size()) {
			const int u = que.front();
			que.pop();
			for(auto i : g[u]) {
				auto [v, c] = e[i];
				if(h[v] == -1 && c > 0) {
					h[v] = h[u] + 1;
					if(v == t) {
						return true;
					}
					que.push(v);
				}
			}
		}
		return false;
	}

	T dfs(int u, int t, T f) {
		if(u == t) {
			return f;
		}
		auto r = f;
		for(int &i = cur[u]; i < (int)g[u].size(); i++) {
			const int j = g[u][i];
			auto [v, c] = e[j];
			if(c > 0 && h[v] == h[u] + 1) {
				auto a = dfs(v, t, std::min(r, c));
				e[j].cap -= a;
				e[j ^ 1].cap += a;
				r -= a;
				if(r == 0) {
					return f;
				}
			}
		}
		return f - r;
	}

	void addadge(int u, int v, T c) {
		g[u].push_back(e.size());
		e.emplace_back(v, c);
		g[v].push_back(e.size());
		e.emplace_back(u, 0); 
	} // 反向边用 j ^ 1 求

	T flow(int s, int t) {
		T ans = 0;
		while(bfs(s, t)) {
			cur.assign(n + 1, 0);
			ans += dfs(s, t, std::numeric_limits<T>::max());
		}
		return ans;
	}

	std::vector<int> minCut() {
		std::vector<int> a(n + 1);
		for(int i = 1; i <= n; i++) {
			a[i] = (h[i] != -1);
		}
		return a;
	}
};

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	std::cout.tie(nullptr);

	int m, n;
	std::string ss;
	std::cin >> m >> n;
	std::getline(std::cin, ss);
	int s = n + m + 1, t = n + m + 2;
	Maxflow<i64> G(t + 10);
	i64 sum = 0;
	for(int i = 1; i <= m; i++) {
		std::getline(std::cin, ss);
		std::stringstream sss;
		sss << ss; 
		i64 x;
		sss >> x;
		sum += x;
		G.addadge(s, i, x);
		while(!sss.eof()) {
			sss >> x;
			G.addadge(i, x + m, 1e18);
		}
	}
	for(int i = 1; i <= n; i++) {
		i64 x;
		std::cin >> x;
		G.addadge(i + m, t, x);
	}
	auto mc = G.flow(s, t);
	auto ans = G.minCut();
	for(int i = 1; i <= m; i++) {
		if(ans[i]) {
			std::cout << i << ' ';
		}
	} std::cout << '\n';
	for(int i = m + 1; i <= m + n; i++) {
		if(ans[i]) {
			std::cout << i - m << ' ';
		}
	} std::cout << '\n';
	std::cout << sum - mc << '\n';
}
```

## mint/comb

```cpp
using i64 = long long;

template<class T>
constexpr T power(T a, i64 b) {
    T res = 1;
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
}
 
constexpr i64 mul(i64 a, i64 b, i64 p) {
    i64 res = a * b - i64(1.L * a * b / p) * p;
    res %= p;
    if (res < 0) {
        res += p;
    }
    return res;
}
template<i64 P>
struct MLong {
    i64 x;
    constexpr MLong() : x{} {}
    constexpr MLong(i64 x) : x{norm(x % getMod())} {}
    
    static i64 Mod;
    constexpr static i64 getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(i64 Mod_) {
        Mod = Mod_;
    }
    constexpr i64 norm(i64 x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr i64 val() const {
        return x;
    }
    explicit constexpr operator i64() const {
        return x;
    }
    constexpr MLong operator-() const {
        MLong res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MLong inv() const {
        assert(x != 0);
        return power(*this, getMod() - 2);
    }
    constexpr MLong &operator*=(MLong rhs) & {
        x = mul(x, rhs.x, getMod());
        return *this;
    }
    constexpr MLong &operator+=(MLong rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MLong &operator-=(MLong rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MLong &operator/=(MLong rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MLong operator*(MLong lhs, MLong rhs) {
        MLong res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MLong operator+(MLong lhs, MLong rhs) {
        MLong res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MLong operator-(MLong lhs, MLong rhs) {
        MLong res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MLong operator/(MLong lhs, MLong rhs) {
        MLong res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MLong &a) {
        i64 v;
        is >> v;
        a = MLong(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MLong &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MLong lhs, MLong rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MLong lhs, MLong rhs) {
        return lhs.val() != rhs.val();
    }
};

template<>
i64 MLong<0LL>::Mod = i64(1E18) + 9;

template<int P>
struct MInt {
    int x;
    constexpr MInt() : x{} {}
    constexpr MInt(i64 x) : x{norm(x % getMod())} {}
    
    static int Mod;
    constexpr static int getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(int Mod_) {
        Mod = Mod_;
    }
    constexpr int norm(int x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr int val() const {
        return x;
    }
    explicit constexpr operator int() const {
        return x;
    }
    constexpr MInt operator-() const {
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv() const {
        assert(x != 0);
        return power(*this, getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs) & {
        x = 1LL * x * rhs.x % getMod();
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MInt &a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MInt &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs, MInt rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs, MInt rhs) {
        return lhs.val() != rhs.val();
    }
};

template<>
int MInt<0>::Mod = 998244353;

template<int V, int P>
constexpr MInt<P> CInv = MInt<P>(V).inv();

constexpr int P = 998244353;
using Z = MInt<P>;

struct Comb {
    int n;
    std::vector<Z> _fac;
    std::vector<Z> _invfac;
    std::vector<Z> _inv;
    
    Comb() : n{0}, _fac{1}, _invfac{1}, _inv{0} {}
    Comb(int n) : Comb() {
        init(n);
    }
    
    void init(int m) {
        m = std::min(m, Z::getMod() - 1);
        if (m <= n) return;
        _fac.resize(m + 1);
        _invfac.resize(m + 1);
        _inv.resize(m + 1);
        
        for (int i = n + 1; i <= m; i++) {
            _fac[i] = _fac[i - 1] * i;
        }
        _invfac[m] = _fac[m].inv();
        for (int i = m; i > n; i--) {
            _invfac[i - 1] = _invfac[i] * i;
            _inv[i] = _invfac[i] * _fac[i - 1];
        }
        n = m;
    }
    
    Z fac(int m) {
        if (m > n) init(2 * m);
        return _fac[m];
    }
    Z invfac(int m) {
        if (m > n) init(2 * m);
        return _invfac[m];
    }
    Z inv(int m) {
        if (m > n) init(2 * m);
        return _inv[m];
    }
    Z binom(int n, int m) {
        if (n < m || m < 0) return 0;
        return fac(n) * invfac(m) * invfac(n - m);
    }
} comb;

```

# 三、计算几何

## cal_geo

```cpp
template<class T>
struct Point {
    T x;
    T y;
    Point(T x_ = 0, T y_ = 0) : x(x_), y(y_) {}
     
    template<class U>
    operator Point<U>() {
        return Point<U>(U(x), U(y));
    }
    Point &operator+=(Point p) & {
        x += p.x;
        y += p.y;
        return *this;
    }
    Point &operator-=(Point p) & {
        x -= p.x;
        y -= p.y;
        return *this;
    }
    Point &operator*=(T v) & {
        x *= v;
        y *= v;
        return *this;
    }
    Point operator-() const {
        return Point(-x, -y);
    }
    friend Point operator+(Point a, Point b) {
        return a += b;
    }
    friend Point operator-(Point a, Point b) {
        return a -= b;
    }
    friend Point operator*(Point a, T b) {
        return a *= b;
    }
    friend Point operator*(T a, Point b) {
        return b *= a;
    }
    friend bool operator==(Point a, Point b) {
        return a.x == b.x && a.y == b.y;
    }
    friend std::istream &operator>>(std::istream &is, Point &p) {
        return is >> p.x >> p.y;
    }
    friend std::ostream &operator<<(std::ostream &os, Point p) {
        return os << "(" << p.x << ", " << p.y << ")";
    }
};
 
template<class T>
T dot(Point<T> a, Point<T> b) {
    return a.x * b.x + a.y * b.y;
}
 
template<class T>
T cross(Point<T> a, Point<T> b) {
    return a.x * b.y - a.y * b.x;
}
 
template<class T>
T square(Point<T> p) {
    return dot(p, p);
}
 
template<class T>
double length(Point<T> p) {
    return std::sqrt(double(square(p)));
}
 
long double length(Point<long double> p) {
    return std::sqrt(square(p));
}
 
template<class T>
struct Line {
    Point<T> a;
    Point<T> b;
    Line(Point<T> a_ = Point<T>(), Point<T> b_ = Point<T>()) : a(a_), b(b_) {}
};
 
template<class T>
Point<T> rotate(Point<T> a) {
    return Point(-a.y, a.x);
}
 
template<class T>
int sgn(Point<T> a) {
    return a.y > 0 || (a.y == 0 && a.x > 0) ? 1 : -1;
}
 
template<class T>
bool pointOnLineLeft(Point<T> p, Line<T> l) {
    return cross(l.b - l.a, p - l.a) > 0;
}
 
template<class T>
Point<T> lineIntersection(Line<T> l1, Line<T> l2) {
    return l1.a + (l1.b - l1.a) * (cross(l2.b - l2.a, l1.a - l2.a) / cross(l2.b - l2.a, l1.a - l1.b));
}
 
template<class T>
bool pointOnSegment(Point<T> p, Line<T> l) {
    return cross(p - l.a, l.b - l.a) == 0 && std::min(l.a.x, l.b.x) <= p.x && p.x <= std::max(l.a.x, l.b.x)
        && std::min(l.a.y, l.b.y) <= p.y && p.y <= std::max(l.a.y, l.b.y);
}
 
template<class T>
bool pointInPolygon(Point<T> a, std::vector<Point<T>> p) {
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (pointOnSegment(a, Line(p[i], p[(i + 1) % n]))) {
            return true;
        }
    }
     
    int t = 0;
    for (int i = 0; i < n; i++) {
        auto u = p[i];
        auto v = p[(i + 1) % n];
        if (u.x < a.x && v.x >= a.x && pointOnLineLeft(a, Line(v, u))) {
            t ^= 1;
        }
        if (u.x >= a.x && v.x < a.x && pointOnLineLeft(a, Line(u, v))) {
            t ^= 1;
        }
    }
     
    return t == 1;
}
 
// 0 : not intersect
// 1 : strictly intersect
// 2 : overlap
// 3 : intersect at endpoint
template<class T>
std::tuple<int, Point<T>, Point<T>> segmentIntersection(Line<T> l1, Line<T> l2) {
    if (std::max(l1.a.x, l1.b.x) < std::min(l2.a.x, l2.b.x)) {
        return {0, Point<T>(), Point<T>()};
    }
    if (std::min(l1.a.x, l1.b.x) > std::max(l2.a.x, l2.b.x)) {
        return {0, Point<T>(), Point<T>()};
    }
    if (std::max(l1.a.y, l1.b.y) < std::min(l2.a.y, l2.b.y)) {
        return {0, Point<T>(), Point<T>()};
    }
    if (std::min(l1.a.y, l1.b.y) > std::max(l2.a.y, l2.b.y)) {
        return {0, Point<T>(), Point<T>()};
    }
    if (cross(l1.b - l1.a, l2.b - l2.a) == 0) {
        if (cross(l1.b - l1.a, l2.a - l1.a) != 0) {
            return {0, Point<T>(), Point<T>()};
        } else {
            auto maxx1 = std::max(l1.a.x, l1.b.x);
            auto minx1 = std::min(l1.a.x, l1.b.x);
            auto maxy1 = std::max(l1.a.y, l1.b.y);
            auto miny1 = std::min(l1.a.y, l1.b.y);
            auto maxx2 = std::max(l2.a.x, l2.b.x);
            auto minx2 = std::min(l2.a.x, l2.b.x);
            auto maxy2 = std::max(l2.a.y, l2.b.y);
            auto miny2 = std::min(l2.a.y, l2.b.y);
            Point<T> p1(std::max(minx1, minx2), std::max(miny1, miny2));
            Point<T> p2(std::min(maxx1, maxx2), std::min(maxy1, maxy2));
            if (!pointOnSegment(p1, l1)) {
                std::swap(p1.y, p2.y);
            }
            if (p1 == p2) {
                return {3, p1, p2};
            } else {
                return {2, p1, p2};
            }
        }
    }
    auto cp1 = cross(l2.a - l1.a, l2.b - l1.a);
    auto cp2 = cross(l2.a - l1.b, l2.b - l1.b);
    auto cp3 = cross(l1.a - l2.a, l1.b - l2.a);
    auto cp4 = cross(l1.a - l2.b, l1.b - l2.b);
     
    if ((cp1 > 0 && cp2 > 0) || (cp1 < 0 && cp2 < 0) || (cp3 > 0 && cp4 > 0) || (cp3 < 0 && cp4 < 0)) {
        return {0, Point<T>(), Point<T>()};
    }
     
    Point p = lineIntersection(l1, l2);
    if (cp1 != 0 && cp2 != 0 && cp3 != 0 && cp4 != 0) {
        return {1, p, p};
    } else {
        return {3, p, p};
    }
}
 
template<class T>
bool segmentInPolygon(Line<T> l, std::vector<Point<T>> p) {
    int n = p.size();
    if (!pointInPolygon(l.a, p)) {
        return false;
    }
    if (!pointInPolygon(l.b, p)) {
        return false;
    }
    for (int i = 0; i < n; i++) {
        auto u = p[i];
        auto v = p[(i + 1) % n];
        auto w = p[(i + 2) % n];
        auto [t, p1, p2] = segmentIntersection(l, Line(u, v));
         
        if (t == 1) {
            return false;
        }
        if (t == 0) {
            continue;
        }
        if (t == 2) {
            if (pointOnSegment(v, l) && v != l.a && v != l.b) {
                if (cross(v - u, w - v) > 0) {
                    return false;
                }
            }
        } else {
            if (p1 != u && p1 != v) {
                if (pointOnLineLeft(l.a, Line(v, u))
                    || pointOnLineLeft(l.b, Line(v, u))) {
                    return false;
                }
            } else if (p1 == v) {
                if (l.a == v) {
                    if (pointOnLineLeft(u, l)) {
                        if (pointOnLineLeft(w, l)
                            && pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    } else {
                        if (pointOnLineLeft(w, l)
                            || pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    }
                } else if (l.b == v) {
                    if (pointOnLineLeft(u, Line(l.b, l.a))) {
                        if (pointOnLineLeft(w, Line(l.b, l.a))
                            && pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    } else {
                        if (pointOnLineLeft(w, Line(l.b, l.a))
                            || pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    }
                } else {
                    if (pointOnLineLeft(u, l)) {
                        if (pointOnLineLeft(w, Line(l.b, l.a))
                            || pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    } else {
                        if (pointOnLineLeft(w, l)
                            || pointOnLineLeft(w, Line(u, v))) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}
 
template<class T>
std::vector<Point<T>> hp(std::vector<Line<T>> lines) {
    std::sort(lines.begin(), lines.end(), [&](auto l1, auto l2) {
        auto d1 = l1.b - l1.a;
        auto d2 = l2.b - l2.a;
         
        if (sgn(d1) != sgn(d2)) {
            return sgn(d1) == 1;
        }
         
        return cross(d1, d2) > 0;
    });
     
    std::deque<Line<T>> ls;
    std::deque<Point<T>> ps;
    for (auto l : lines) {
        if (ls.empty()) {
            ls.push_back(l);
            continue;
        }
         
        while (!ps.empty() && !pointOnLineLeft(ps.back(), l)) {
            ps.pop_back();
            ls.pop_back();
        }
         
        while (!ps.empty() && !pointOnLineLeft(ps[0], l)) {
            ps.pop_front();
            ls.pop_front();
        }
         
        if (cross(l.b - l.a, ls.back().b - ls.back().a) == 0) {
            if (dot(l.b - l.a, ls.back().b - ls.back().a) > 0) {
                 
                if (!pointOnLineLeft(ls.back().a, l)) {
                    assert(ls.size() == 1);
                    ls[0] = l;
                }
                continue;
            }
            return {};
        }
         
        ps.push_back(lineIntersection(ls.back(), l));
        ls.push_back(l);
    }
     
    while (!ps.empty() && !pointOnLineLeft(ps.back(), ls[0])) {
        ps.pop_back();
        ls.pop_back();
    }
    if (ls.size() <= 2) {
        return {};
    }
    ps.push_back(lineIntersection(ls[0], ls.back()));
     
    return std::vector(ps.begin(), ps.end());
}
```

## cal_geo_du

```cpp
const double EPS = 1e-9;
//由于硬件限制，浮点数运算有误差，eps用来消除误差
inline int sign(double a) { return a < -EPS ? -1 : a > EPS; }
//判断数符号，负数返回-1，0返回0，正数返回1
inline int cmp(double a, double b) { return sign(a - b); }
//比较两数大小
//点类，向量类
//因为有许多操作相似，所以并在一起
struct P {
    double x, y;
    //点表示坐标，向量表示向量
    P() {}
    P(double _x, double _y) : x(_x), y(_y) {}
    //构造函数
    P operator+(P p) { return {x + p.x, y + p.y}; }
    P operator-(P p) { return {x - p.x, y - p.y}; }
    P operator*(double d) { return {x * d, y * d}; }
    P operator/(double d) { return {x / d, y / d}; }
    //向量加减乘除
    bool operator<(P p) const
    {
        int c = cmp(x, p.x);
        if (c)
            return c == -1;
        return cmp(y, p.y) == -1;
    }
    bool operator==(P o) const
    {
        return cmp(x, o.x) == 0 && cmp(y, o.y) == 0;
    }
    //比较字典序
    double dot(P p) { return x * p.x + y * p.y; }
    //点积
    double det(P p) { return x * p.y - y * p.x; }
    //叉积
    double distTo(P p) { return (*this - p).abs(); }
    //点距离
    double alpha() { return atan2(y, x); }
    void read() { cin >> x >> y; }
    void write() { cout << "(" << x << "," << y << ")" << endl; }
    double abs() { return sqrt(abs2()); }
    double abs2() { return x * x + y * y; }
    P rot90() { return P(-y, x); }
    P unit() { return *this / abs(); }
    int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }
    //判断点在极角坐标系上半边还是下半边，极点和极轴也算上半边
    P rot(double an) { return {x * cos(an) - y * sin(an), x * sin(an) + y * cos(an)}; }
    //向量旋转
};
//线类，半平面类
struct L { // ps[0] -> ps[1]
    P ps[2];
    P &operator[](int i) { return ps[i]; }
    P dir() { return ps[1] - ps[0]; }
    L(P a, P b)
    {
        ps[0] = a;
        ps[1] = b;
    }
    bool include(P p) { return sign((ps[1] - ps[0]).det(p - ps[0])) > 0; }
    L push()
    { // push eps outward
        const double eps = 1e-8;
        P delta = (ps[1] - ps[0]).rot90().unit() * eps;
        return {ps[0] + delta, ps[1] + delta};
    }
};

#define cross(p1, p2, p3) ((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y))
#define crossOp(p1, p2, p3) sign(cross(p1, p2, p3))
//叉积
bool chkLL(P p1, P p2, P q1, P q2) {
    double a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return sign(a1 + a2) != 0;
}
//判断向量平行
P isLL(P p1, P p2, P q1, P q2) {
    double a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return (p1 * a2 + p2 * a1) / (a1 + a2);
}
P isLL(L l1, L l2) { return isLL(l1[0], l1[1], l2[0], l2[1]); }
//求直线交点
bool intersect(double l1, double r1, double l2, double r2) {
    if (l1 > r1)
        swap(l1, r1);
    if (l2 > r2)
        swap(l2, r2);
    return !(cmp(r1, l2) == -1 || cmp(r2, l1) == -1);
}
bool isSS(P p1, P p2, P q1, P q2) {
    return intersect(p1.x, p2.x, q1.x, q2.x) && intersect(p1.y, p2.y, q1.y, q2.y) &&
           crossOp(p1, p2, q1) * crossOp(p1, p2, q2) <= 0 && crossOp(q1, q2, p1) * crossOp(q1, q2, p2) <= 0;
}

bool isSS_strict(P p1, P p2, P q1, P q2) {
    return crossOp(p1, p2, q1) * crossOp(p1, p2, q2) < 0 && crossOp(q1, q2, p1) * crossOp(q1, q2, p2) < 0;
}
//判断线段相交，交在端点算不算分为严格不严格
bool isMiddle(double a, double m, double b)
{
    return sign(a - m) == 0 || sign(b - m) == 0 || (a < m != b < m);
}
bool isMiddle(P a, P m, P b)
{
    return isMiddle(a.x, m.x, b.x) && isMiddle(a.y, m.y, b.y);
}
bool onSeg(P p1, P p2, P q)
{
    return crossOp(p1, p2, q) == 0 && isMiddle(p1, q, p2);
}
bool onSeg_strict(P p1, P p2, P q)
{
    return crossOp(p1, p2, q) == 0 && sign((q - p1).dot(p1 - p2)) * sign((q - p2).dot(p1 - p2)) < 0;
}
//点在线段上判定
P proj(P p1, P p2, P q)
{
    P dir = p2 - p1;
    return p1 + dir * (dir.dot(q - p1) / dir.abs2());
}
P reflect(P p1, P p2, P q)
{
    return proj(p1, p2, q) * 2 - q;
}
double nearest(P p1, P p2, P q)
{
    if (p1 == p2)
        return p1.distTo(q);
    P h = proj(p1, p2, q);
    if (isMiddle(p1, h, p2))
        return q.distTo(h);
    return min(p1.distTo(q), p2.distTo(q));
}
//投影，反射，最近点
//最近点是线段外一点到线段上的点的最短距离
double disSS(P p1, P p2, P q1, P q2)
{
    if (isSS(p1, p2, q1, q2))
        return 0;
    return min(min(nearest(p1, p2, q1), nearest(p1, p2, q2)), min(nearest(q1, q2, p1), nearest(q1, q2, p2)));
}
//线段距离
double rad(P p1, P p2)
{
    return atan2l(p1.det(p2), p1.dot(p2));
}

double incircle(P p1, P p2, P p3)
{
    double A = p1.distTo(p2);
    double B = p2.distTo(p3);
    double C = p3.distTo(p1);
    return sqrtl(A * B * C / (A + B + C));
}

// polygon
//简单多边形的问题只有判断点在多边形内，和多边形面积简单，其他只做凸多边形
double area(vector<P> ps)
{
    double ret = 0;
    for(int i=0;i<ps.size();++i) 
        ret += ps[i].det(ps[(i + 1) % ps.size()]);
    return ret / 2;
}
//多边形面积
int contain(vector<P> ps, P p)
{ // 2:inside,1:on_seg,0:outside
    int n = ps.size(), ret = 0;
    for(int i = 0; i < n; i++)
    {
        P u = ps[i], v = ps[(i + 1) % n];
        if (onSeg(u, v, p))
            return 1;
        if (cmp(u.y, v.y) <= 0)
            swap(u, v);
        if (cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0)
            continue;
        ret ^= crossOp(p, u, v) > 0;
    }
    return ret * 2;
}
//判断点在多边形内
vector<P> convexHull(vector<P> ps)
{
    int n = ps.size();
    if (n <= 1)
        return ps;
    sort(ps.begin(), ps.end());
    vector<P> qs(n * 2);
    int k = 0;
    for (int i = 0; i < n; qs[k++] = ps[i++])
        while (k > 1 && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0)
            --k;
    for (int i = n - 2, t = k; i >= 0; qs[k++] = ps[i--])
        while (k > t && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0)
            --k;
    qs.resize(k - 1);
    return qs;
}
vector<P> convexHullNonStrict(vector<P> ps)
{
    // caution: need to unique the Ps first
    map<P, int> mp;
    for(auto x : ps) {
    	mp[x]++;
    }
    vector<P> pp, pq;
    for(auto [x, y] : mp) {
    	pp.push_back(x);
    }
    int n = pp.size();
    if (n <= 1)
        return pp;
    sort(pp.begin(), pp.end());
    vector<P> qs(n * 2);
    int k = 0;
    for (int i = 0; i < n; qs[k++] = pp[i++])
        while (k > 1 && crossOp(qs[k - 2], qs[k - 1], pp[i]) < 0)
            --k;
    for (int i = n - 2, t = k; i >= 0; qs[k++] = pp[i--])
        while (k > t && crossOp(qs[k - 2], qs[k - 1], pp[i]) < 0)
            --k;
    qs.resize(k - 1);
 
    for(auto x : qs) {
    	for(int i = 1; i <= mp[x]; i++) {
    		pq.push_back(x);
    	}
    }
    return pq;
}
//凸包

double convexDiameter(vector<P> ps)
{
    int n = ps.size();
    if (n <= 1)
        return 0;
    int is = 0, js = 0;
    for(int k = 1; k < n; k++) is = ps[k] < ps[is] ? k : is, js = ps[js] < ps[k] ? k : js;
    int i = is, j = js;
    double ret = ps[i].distTo(ps[j]);
    do
    {
        if ((ps[(i + 1) % n] - ps[i]).det(ps[(j + 1) % n] - ps[j]) >= 0)
            (++j) %= n;
        else
            (++i) %= n;
        ret = max(ret, ps[i].distTo(ps[j]));
    } while (i != is || j != js);
    return ret;
}
//凸包直径
vector<P> convexCut(const vector<P> &ps, P q1, P q2)
{
    vector<P> qs;
    int n = ps.size();
    for(int i = 0; i < n; i++)
    {
        P p1 = ps[i], p2 = ps[(i + 1) % n];
        int d1 = crossOp(q1, q2, p1), d2 = crossOp(q1, q2, p2);
        if (d1 >= 0)
            qs.push_back(p1);
        if (d1 * d2 < 0)
            qs.push_back(isLL(p1, p2, q1, q2));
    }
    return qs;
}
//直线切割凸包，返回直线左边凸包的点
double min_dist(vector<P> &ps, int l, int r)
{
    if (r - l <= 5)
    {
        double ret = 1e18;
        for(int i=l;i<r;++i)
            for(int j=l;j<i;++j)
                ret = min(ret, ps[i].distTo(ps[j]));
        return ret;
    }
    int m = (l + r) >> 1;
    double ret = min(min_dist(ps, l, m), min_dist(ps, m, r));
    vector<P> qs;
    for(int i=l;i<r;++i)
        if (abs(ps[i].x - ps[m].x) <= ret)
            qs.push_back(ps[i]);
    sort(qs.begin(), qs.end(), [](P a, P b) -> bool
         { return a.y < b.y; });
    for(int i=1;i<qs.size();++i) 
        for (int j = i - 1; j >= 0 && qs[j].y >= qs[i].y - ret; --j)
            ret = min(ret, qs[i].distTo(qs[j]));
    return ret;
}
//平面最近点对,[l,r)，要求ps按x升序
int type(P o1, double r1, P o2, double r2)
{
    double d = o1.distTo(o2);
    if (cmp(d, r1 + r2) == 1)
        return 4;
    if (cmp(d, r1 + r2) == 0)
        return 3;
    if (cmp(d, abs(r1 - r2)) == 1)
        return 2;
    if (cmp(d, abs(r1 - r2)) == 0)
        return 1;
    return 0;
}

vector<P> isCL(P o, double r, P p1, P p2)
{
    if (cmp(abs((o - p1).det(p2 - p1) / p1.distTo(p2)), r) > 0)
        return {};
    double x = (p1 - o).dot(p2 - p1), y = (p2 - p1).abs2(), d = x * x - y * ((p1 - o).abs2() - r * r);
    d = max(d, (double)0.0);
    P m = p1 - (p2 - p1) * (x / y), dr = (p2 - p1) * (sqrt(d) / y);
    return {m - dr, m + dr}; // along dir: p1->p2
}

vector<P> isCC(P o1, double r1, P o2, double r2)
{ // need to check whether two circles are the same
    double d = o1.distTo(o2);
    if (cmp(d, r1 + r2) == 1)
        return {};
    if (cmp(d, abs(r1 - r2)) == -1)
        return {};
    d = min(d, r1 + r2);
    double y = (r1 * r1 + d * d - r2 * r2) / (2 * d), x = sqrt(r1 * r1 - y * y);
    P dr = (o2 - o1).unit();
    P q1 = o1 + dr * y, q2 = dr.rot90() * x;
    return {q1 - q2, q1 + q2}; // along circle 1
}

vector<P> tanCP(P o, double r, P p)
{
    double x = (p - o).abs2(), d = x - r * r;
    if (sign(d) <= 0)
        return {}; // on circle => no tangent
    P q1 = o + (p - o) * (r * r / x);
    P q2 = (p - o).rot90() * (r * sqrt(d) / x);
    return {q1 - q2, q1 + q2}; // counter clock-wise
}

vector<L> extanCC(P o1, double r1, P o2, double r2)
{
    vector<L> ret;
    if (cmp(r1, r2) == 0)
    {
        P dr = (o2 - o1).unit().rot90() * r1;
        ret.push_back(L(o1 + dr, o2 + dr)), ret.push_back(L(o1 - dr, o2 - dr));
    }
    else
    {
        P p = (o2 * r1 - o1 * r2) / (r1 - r2);
        vector<P> ps = tanCP(o1, r1, p), qs = tanCP(o2, r2, p);
        for(int i = 0; i < min(ps.size(), qs.size()); i++) ret.push_back(L(ps[i], qs[i])); // c1 counter-clock wise
    }
    return ret;
}

vector<L> intanCC(P o1, double r1, P o2, double r2)
{
    vector<L> ret;
    P p = (o1 * r2 + o2 * r1) / (r1 + r2);
    vector<P> ps = tanCP(o1, r1, p), qs = tanCP(o2, r2, p);
    for(int i = 0; i < min(ps.size(), qs.size()); i++) ret.push_back(L(ps[i], qs[i])); // c1 counter-clock wise
    return ret;
}

double areaCT(double r, P p1, P p2)
{
    vector<P> is = isCL(P(0, 0), r, p1, p2);
    if (is.empty())
        return r * r * rad(p1, p2) / 2;
    bool b1 = cmp(p1.abs2(), r * r) == 1, b2 = cmp(p2.abs2(), r * r) == 1;
    if (b1 && b2)
    {
        if (sign((p1 - is[0]).dot(p2 - is[0])) <= 0 &&
            sign((p1 - is[0]).dot(p2 - is[0])) <= 0)
            return r * r * (rad(p1, is[0]) + rad(is[1], p2)) / 2 + is[0].det(is[1]) / 2;
        else
            return r * r * rad(p1, p2) / 2;
    }
    if (b1)
        return (r * r * rad(p1, is[0]) + is[0].det(p2)) / 2;
    if (b2)
        return (p1.det(is[1]) + r * r * rad(is[1], p2)) / 2;
    return p1.det(p2) / 2;
}

bool parallel(L l0, L l1) { return sign(l0.dir().det(l1.dir())) == 0; }
bool cmp(P a, P b)
{
    if (a.quad() != b.quad())
    {
        return a.quad() < b.quad();
    }
    else
    {
        return sign(a.det(b)) > 0;
    }
}
//极角排序
bool sameDir(L l0, L l1) { return parallel(l0, l1) && sign(l0.dir().dot(l1.dir())) == 1; }
bool operator<(L l0, L l1)
{
    if (sameDir(l0, l1))
    {
        return l1.include(l0[0]);
    }
    else
    {
        return cmp(l0.dir(), l1.dir());
    }
}
bool check(L u, L v, L w)
{
    return w.include(isLL(u, v));
}
vector<P> halfPlaneIS(vector<L> &l)
{
    sort(l.begin(), l.end());
    deque<L> q;
    for (int i = 0; i < (int)l.size(); ++i)
    {
        if (i && sameDir(l[i], l[i - 1]))
            continue;
        while (q.size() > 1 && !check(q[q.size() - 2], q[q.size() - 1], l[i]))
            q.pop_back();
        while (q.size() > 1 && !check(q[1], q[0], l[i]))
            q.pop_front();
        q.push_back(l[i]);
    }
    while (q.size() > 2 && !check(q[q.size() - 2], q[q.size() - 1], q[0]))
        q.pop_back();
    while (q.size() > 2 && !check(q[1], q[0], q[q.size() - 1]))
        q.pop_front();
    vector<P> ret;
    for (int i = 0; i < (int)q.size(); ++i)
        ret.push_back(isLL(q[i], q[(i + 1) % q.size()]));
    return ret;
}
//半平面交
P inCenter(P A, P B, P C)
{
    double a = (B - C).abs(), b = (C - A).abs(), c = (A - B).abs();
    return (A * a + B * b + C * c) / (a + b + c);
}
//内心，角平分线的交点
P circumCenter(P a, P b, P c)
{
    P bb = b - a, cc = c - a;
    double db = bb.abs2(), dc = cc.abs2(), d = 2 * bb.det(cc);
    return a - P(bb.y * dc - cc.y * db, cc.x * db - bb.x * dc) / d;
}
//外心，垂直平分线的交点
P orthoCenter(P a, P b, P c)
{
    P ba = b - a, ca = c - a, bc = b - c;
    double Y = ba.y * ca.y * bc.y,
           A = ca.x * ba.y - ba.x * ca.y,
           x0 = (Y + ca.x * ba.y * b.x - ba.x * ca.y * c.x) / A,
           y0 = -ba.x * (x0 - c.x) / ba.y + ca.y;
    return {x0, y0};
}
//垂心，垂线的交点
```

# 四、数据结构

## 线段树加乘

```cpp
#include<bits/stdc++.h>
using namespace std;
#define int long long
#define AC return
#define Please 0

const int N = 1e5 + 10;
int n, m, p;
int x, y, z, k;

struct Node {
	int l, r;
	int lazy1, lazy2;
		// 1 cheng 2 jia
	int sum, sz;
	void fuc() {
		lazy1 %= p;
		lazy2 %= p;
		sum %= p;
	}
}tr[N * 4];

void pushup(int u) {
	tr[u].sum = (tr[u << 1].sum + tr[u << 1 | 1].sum) % p;
	tr[u].fuc();
}

void pushdown(int u) {
	auto &root = tr[u], &left = tr[u << 1], &right = tr[u << 1| 1];
	left.sum *= root.lazy1, right.sum *= root.lazy1;
	left.lazy1 *= root.lazy1, right.lazy1 *= root.lazy1;
	left.lazy2 *= root.lazy1, right.lazy2 *= root.lazy1;
	root.lazy1 = 1;
	left.sum += root.lazy2 * left.sz, right.sum += root.lazy2 * right.sz;
	left.lazy2 += root.lazy2, right.lazy2 += root.lazy2;
	root.lazy2 = 0;
	root.fuc(), left.fuc(), right.fuc();
}

void build(int u, int l, int r) {
	auto &root = tr[u];
	root.l = l, root.r = r;
	root.sz = r - l + 1;
	root.sum = 0;
	root.lazy1 = 1;
	root.lazy2 = 0;
	if(l == r) {
		return;
	}
	int mid = (l + r) >> 1;
	build(u << 1, l, mid);
	build(u << 1 | 1, mid + 1, r);
}

void modify(int u, int l, int r, int x, int op) {
	int L = tr[u].l, R = tr[u].r;
	int mid = (L + R) >> 1;
	if(l <= L && R <= r) {
		if(op == 1) {
			tr[u].lazy1 *= x;
			tr[u].lazy2 *= x;
			tr[u].sum *= x;
			tr[u].fuc();
		} else {
			tr[u].lazy2 += x;
			tr[u].sum += x * tr[u].sz;
			tr[u].fuc();
		}
		return;
	}
	pushdown(u);
	if(l <= mid) {
		modify(u << 1, l, r, x, op);
	}
	if(r > mid) {
		modify(u << 1 | 1, l, r, x, op);
	}
	pushup(u);
}

int query(int u, int l, int r) {
	int L = tr[u].l, R = tr[u].r;
	int ans = 0;
	if(l <= L && R <= r) {
		return tr[u].sum;
	}
	int mid = (L + R) >> 1;
	pushdown(u);
	if(l <= mid) {
		ans += query(u << 1, l, r);
		ans %= p;
	}
	if(r > mid) {
		ans += query(u << 1 | 1, l, r);
		ans %= p;
	}
	return ans;
}

void solve() {
	cin >> n >> m >> p;
	build(1, 1, n);
	for(int i = 1; i <= n; i++) {
		cin >> x;
		modify(1, i, i, x, 2);
	}
	while(m--) {
		int op;
		cin >> op;
		if(op != 3) {
			cin >> x >> y >> k;
			modify(1, x, y, k, op);
            // op == 1 plus
            // op == 2 multiply
		} else {
			cin >> x >> y;
			cout << query(1, x, y) << endl;
		}
	}
}

signed main() {
	solve();
	AC Please;
}
```



## jlyseg_template

```cpp
#include <bits/stdc++.h>
 
using i64 = long long;
template<class Info, class Tag>
struct LazySegmentTree {
    int n;
    std::vector<Info> info;
    std::vector<Tag> tag;
    LazySegmentTree() : n(0) {}
    LazySegmentTree(int n_, Info v_ = Info()) {
        init(n_, v_);
    }
    template<class T>
    LazySegmentTree(std::vector<T> init_) {
        init(init_);
    }
    void init(int n_, Info v_ = Info()) {
        init(std::vector(n_, v_));
    }
    template<class T>
    void init(std::vector<T> init_) {
        n = init_.size();
        info.assign(4 << std::__lg(n), Info());
        tag.assign(4 << std::__lg(n), Tag());
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init_[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void apply(int p, const Tag &v) {
        info[p].apply(v);
        tag[p].apply(v);
    }
    void push(int p) {
        apply(2 * p, tag[p]);
        apply(2 * p + 1, tag[p]);
        tag[p] = Tag();
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        push(p);
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        push(p);
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    void rangeApply(int p, int l, int r, int x, int y, const Tag &v) {
        if (l >= y || r <= x) {
            return;
        }
        if (l >= x && r <= y) {
            apply(p, v);
            return;
        }
        int m = (l + r) / 2;
        push(p);
        rangeApply(2 * p, l, m, x, y, v);
        rangeApply(2 * p + 1, m, r, x, y, v);
        pull(p);
    }
    void rangeApply(int l, int r, const Tag &v) {
        return rangeApply(1, 0, n, l, r, v);
    }
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findLast(2 * p + 1, m, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
    }
};
 
struct Tag {
    i64 a = 0, b = 0;
    void apply(Tag t) {
        a = std::min(a, b + t.a);
        b += t.b;
    }
};
 
int k;
 
struct Info {
    i64 x = 0;
    void apply(Tag t) {
        x += t.a;
        if (x < 0) {
            x = (x % k + k) % k;
        }
        x += t.b - t.a;
    }
};
Info operator+(Info a, Info b) {
    return {a.x + b.x};
}
 
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
     
    int n, m;
    std::cin >> n >> m >> k;
     
    LazySegmentTree<Info, Tag> seg(n);
     
    __int128 ans = 0;
    for (int i = 0; i < m; i++) {
        int o;
        std::cin >> o;
         
        if (o == 1) {
            int l, r, x;
            std::cin >> l >> r >> x;
            l--;
            seg.rangeApply(l, r, {0, x});
            ans += 1LL * (r - l) * x;
        } else {
            int l, r;
            std::cin >> l >> r;
            l--;
            seg.rangeApply(l, r, {-k, -k});
        }
    }
     
    for (int i = 0; i < n; i++) {
        ans -= seg.rangeQuery(i, i + 1).x;
    }
    ans /= k;
     
    std::cout << i64(ans) << "\n";
     
    return 0;
}
```

## jlyFenwick_template

```cpp
#include <bits/stdc++.h>
 
using i64 = long long;
template <typename T>
struct Fenwick {
    int n;
    std::vector<T> a;
     
    Fenwick(int n = 0) {
        init(n);
    }
     
    void init(int n) {
        this->n = n;
        a.assign(n, T());
    }
     
    void add(int x, T v) {
        for (int i = x + 1; i <= n; i += i & -i) {
            a[i - 1] += v;
        }
    }
     
    T sum(int x) {
        auto ans = T();
        for (int i = x; i > 0; i -= i & -i) {
            ans += a[i - 1];
        }
        return ans;
    }
     
    T rangeSum(int l, int r) {
        return sum(r) - sum(l);
    }
     
    int kth(T k) {
        int x = 0;
        for (int i = 1 << std::__lg(n); i; i /= 2) {
            if (x + i <= n && k >= a[x + i - 1]) {
                x += i;
                k -= a[x - 1];
            }
        }
        return x;
    }
};
 
struct Info {
    int x = -1E9;
    Info &operator+=(Info b) & {
        x = std::max(x, b.x);
        return *this;
    }
};
 
 
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
     
    int n;
    std::cin >> n;
     
    std::vector<int> a(n), b(n);
    for (int i = 0; i < n; i++) {
        std::cin >> a[i];
    }
    for (int i = 0; i < n; i++) {
        std::cin >> b[i];
    }
     
    std::vector<int> v;
    for (int i = 0; i < n; i++) {
        v.push_back(std::max(a[i], b[i]));
    }
    std::sort(v.begin(), v.end());
     
    std::vector<std::array<int, 3>> e;
    for (int i = 0; i < n; i++) {
        if (a[i] < b[i]) {
            e.push_back({a[i], b[i], 0});
        }
        if (a[i] > b[i]) {
            e.push_back({b[i], a[i], 1});
        }
    }
     
    std::sort(e.begin(), e.end(), std::greater<std::array<int, 3>>());
     
    int ans = 0;
     
    std::vector fen(2, std::vector(2, Fenwick<Info>(n)));
    for (auto [l, r, t] : e) {
        r = std::lower_bound(v.begin(), v.end(), r) - v.begin();
        ans = std::max(ans, fen[!t][0].sum(r).x);
        ans = std::max(ans, fen[!t][1].sum(n - r).x + v[r]);
        fen[t][0].add(r, Info{v[r] - l});
        fen[t][1].add(n - r - 1, Info{-l});
    }
     
    i64 sum = 0;
    for (int i = 0; i < n; i++) {
        sum += std::abs(a[i] - b[i]);
    }
     
    std::cout << sum - 2LL * ans << "\n";
     
    return 0;
}
```

## 主席树单点维护

```cpp
#include<bits/stdc++.h>

using i64 = long long;

constexpr int N = 2e6 + 10;
struct Node {
	int val;
	int tag;
	int ls, rs;
}tr[N << 4];

int a[N], idx, root[N << 4];

void clone(int u) {
	tr[++idx] = tr[u];
}

void build(int &u, int l, int r) {
	u = ++idx;
	if(l == r) {
		tr[u].val = a[l];
		return;
	} 
	int mid = (l + r) >> 1;
	build(tr[u].ls, l, mid);
	build(tr[u].rs, mid + 1, r);
	return; 
}

void modify(int &u, int v, int l, int r, int beg, int ed, int x) {
	clone(v);
	u = idx;
	if(l == r) {
		tr[u].val = x;
		return;
	}
	int mid = (l + r) >> 1;
	if(beg <= mid) {
		modify(tr[u].ls, tr[v].ls, l, mid, beg, ed, x);
	} else {
		modify(tr[u].rs, tr[v].rs, mid + 1, r, beg, ed, x);
	}
}

void query(int u, int l, int r, int beg, int ed) {
	if(l == r) {
		std::cout << tr[u].val << "\n";
		return;
	}
	int mid = (l + r) >> 1;
	if(beg <= mid) {
		query(tr[u].ls, l, mid, beg, ed);
	} else {
		query(tr[u].rs, mid + 1, r, beg, ed);
	}
}

void solve() {
	int n, m;
	std::cin >> n >> m;
	for(int i = 1; i <= n; i++) {
		std::cin >> a[i];
	}
	build(root[0], 1, n);
	for(int i = 1; i <= m; i++) {
		int v, op, loc, v1;
		std::cin >> v >> op;
		if(op == 1) {
			std::cin >> loc >> v1;
			modify(root[i], root[v], 1, n, loc, loc, v1);
		} else {
			std::cin >> loc;
			query(root[v], 1, n, loc, loc);
			root[i] = root[v];
		}
	}
}

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);

	int t;
	t = 1;
	while(t--) {
		solve();
	}
}
```

<img src="C:\Users\10042\AppData\Roaming\Typora\typora-user-images\image-20231118023848147.png" alt="image-20231118023848147" style="zoom:67%;" />

## 静态区间第 *k* 小

如题，给定 *n* 个整数构成的序列 *a*，将对于指定的闭区间 [*l*, *r*] 查询其区间内的第 *k* 小值。

```cpp
#include<bits/stdc++.h>

using i64 = long long;

constexpr int N = 2e6 + 10;
struct Node {
	int val;
	int tag;
	int ls, rs;
}tr[N << 4];

int a[N], idx, root[N << 4];
std::vector<i64> v;
std::map<int, int> hash, fhash;
void clone(int u) {
	tr[++idx] = tr[u];
}

void pushup(int u, int l, int r) {
	tr[u].val = tr[l].val + tr[r].val;
}

void build(int &u, int l, int r) {
	u = ++idx;
	if(l == r) {
		tr[u].val = 0;
		return;
	} 
	int mid = (l + r) >> 1;
	build(tr[u].ls, l, mid);
	build(tr[u].rs, mid + 1, r);
	pushup(u, tr[u].ls, tr[u].rs);
	return; 
}

void modify(int &u, int v, int l, int r, int beg, int ed) {
	clone(v);
	u = idx;
	if(l == r) {
		tr[u].val++;
		return;
	}
	int mid = (l + r) >> 1;
	if(beg <= mid) {
		modify(tr[u].ls, tr[v].ls, l, mid, beg, ed);
	} else {
		modify(tr[u].rs, tr[v].rs, mid + 1, r, beg, ed);
	}
	pushup(u, tr[u].ls, tr[u].rs);
}

int query(int ql, int qr, int l, int r, int k) {
	if(l == r) {
		return v[l];
	}
	int ze = tr[ql].ls, ye = tr[qr].ls;
	int sum = tr[ye].val - tr[ze].val;
	int mid = (l + r) >> 1;
	if(k <= sum) {
		return query(tr[ql].ls, tr[qr].ls, l, mid, k);
	} else {
		return query(tr[ql].rs, tr[qr].rs, mid + 1, r, k - sum);
	}
}

void solve() {
	int n, m;
	std::cin >> n >> m; 
	v.push_back(-1e18);
	for(int i = 1; i <= n; i++) {
		std::cin >> a[i];
		v.push_back(a[i]);
	}
	std::sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	for(int i = 0; i < v.size(); i++) {
		fhash[v[i]] = i;
	}
	int sz = v.size();
	build(root[0], 1, sz);
	for(int i = 1; i <= n; i++) {
		int v = fhash[a[i]];
		modify(root[i], root[i - 1], 1, sz, v, v);
	}
	for(int i = 1; i <= m; i++) {
		int l, r, k;
		std::cin >> l >> r >> k;
		std::cout << query(root[l - 1], root[r], 1, sz, k) << "\n";
	}
}

int main() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);

	int t;
	t = 1;
	while(t--) {
		solve();
	}
}
```



## splay

![image-20231118024150078](C:\Users\10042\AppData\Roaming\Typora\typora-user-images\image-20231118024150078.png)

```cpp
#include <bits/stdc++.h>
using namespace std;
#define int long long
const int N = 2e5 + 10;

struct Node {
	int s[2];
	int p;
	int v;
	int cnt;
	int sz;
	void init(int p1, int v1) {
		s[0] = s[1] = 0;
		p = p1;
		v = v1;
		cnt = sz = 1;
	}	
} tr[N];

int rt;
int idx;

void pushup(int u) {
	tr[u].sz = tr[tr[u].s[0]].sz + tr[tr[u].s[1]].sz + tr[u].cnt;
}

void rotate(int x) {
	int y = tr[x].p;
	int z = tr[y].p;
	int is = (tr[y].s[1] == x);
	// 
	tr[y].s[is] = tr[x].s[is ^ 1];
	tr[tr[y].s[is]].p = y;
	// 
	tr[x].s[is ^ 1] = y;
	tr[y].p = x;
	// 
	tr[z].s[tr[z].s[1] == y] = x;
	tr[x].p = z;
	// 
	pushup(y);
	pushup(x);
} // 哨兵, -inf, inf

void splay(int x, int k) {
	while(tr[x].p != k) {
		int y = tr[x].p, z = tr[y].p;
		if(z != k) {
			if((tr[y].s[0] == x) ^ (tr[z].s[0] == y)) { // 在同一端
				rotate(x);
			} else {
				rotate(y);
			}
		}
		rotate(x);
	}
	if(k == 0) {
		rt = x;
	}
}

void find(int v) {
	int x = rt;
	while(tr[x].s[v > tr[x].v] && v != tr[x].v) {
		// cout << x << 'x' << tr[x].v << endl;
		x = tr[x].s[v > tr[x].v];
	}
	splay(x, 0);
}

int get_pre(int v) {
	find(v);
	int x = rt;
	if(tr[x].v < v) {
		return x;
	}
	x = tr[x].s[0];
	while(tr[x].s[1]) {
		x = tr[x].s[1];
	}
	return x;
}

int get_suc(int v) {
	find(v);
	// cout << tr[rt].v << endl;
	int x = rt;
	if(tr[x].v > v) {
		return x;
	}
	x = tr[x].s[1];
	while(tr[x].s[0]) {
		x = tr[x].s[0];
	}
	return x;
}

void del(int v) {
	int pre = get_pre(v);
	int suc = get_suc(v);
	splay(pre, 0);
	splay(suc, pre);
	int del = tr[suc].s[0];
	if(tr[del].cnt > 1) {
		tr[del].cnt --;
		splay(del, 0);
	} else if(tr[del].cnt) {
		tr[suc].s[0] = 0;
		splay(suc, 0);
	}
}

int get_rank(int u) {
	find(u);
	return tr[tr[rt].s[0]].sz;
}

int get_val(int k) {
	int x = rt;
	k++;
	while(1) {
		int y = tr[x].s[0];
		if(tr[y].sz + tr[x].cnt < k) {
			k -= tr[y].sz + tr[x].cnt;
			x = tr[x].s[1];
		} else {
			if(tr[y].sz >= k) {
				x = tr[x].s[0];
			} else {
				break;
			}
		}
	}
	splay(x, 0); // 要求
	return tr[x].v;
}

void insert(int jb) {
	int x = rt, p = 0; // p 对应当前节点的父亲
	while(x && tr[x].v != jb) {
		p = x;
		x = tr[x].s[jb > tr[x].v];
	}
	if(x) {
		tr[x].cnt++;
	} else {
		x = ++idx;
		tr[p].s[jb > tr[p].v] = x;
		tr[x].init(p, jb); 
	}
	splay(x, 0);
} 


signed main() {
	int n, x, y;
	cin >> n;
	insert(1e18);
	insert(-1e18);
	while(n--) {
		cin >> x >> y;
		if(x == 1) {
			insert(y);
		}
		if(x == 2) {
			del(y);
		}
		if(x == 3) {
			cout << get_rank(y);
			cout << endl;
		}
		if(x == 4) {
			cout << get_val(y);
			cout << endl;
		}
		if(x == 5) {
			cout << tr[get_pre(y)].v;
			cout << endl; 
		}
		if(x == 6) {
			cout << tr[get_suc(y)].v;
			cout << endl;
		}
	}
}
```



# 五、字符串

## AC自动机

![image-20231118024333907](C:\Users\10042\AppData\Roaming\Typora\typora-user-images\image-20231118024333907.png)

```cpp
#include <bits/stdc++.h>

using namespace std;

struct ACMachine {
	int size;
	struct Node {
		int fail, cnt;
		int next[26]; 
		Node() {
			for(int i = 0; i < 26; i++) {
				next[i] = 0;
			}
			fail = 0;
			cnt = 0;
		}
	};
	vector<Node> tr;
	queue<int> q;
	void init() {
		tr.push_back(Node());
		size = 0;
	}
	void build(string s) {
		int len = s.size(), u = 0;
		s = ' ' + s;
		for(int i = 1; i <= len; i++) {
			int to = s[i] - 'a';
			if(!tr[u].next[to]) {
				tr.push_back(Node());
				size++;
				tr[u].next[to] = size;
			}
			u = tr[u].next[to];
		}
		tr[u].cnt++;
	}
	void getFail() {
		tr[0].fail = 0;
		for(int i = 0; i < 26; i++) {
			if(tr[0].next[i] != 0) {
				tr[tr[0].next[i]].fail = 0;
				q.push(tr[0].next[i]);
			}
		}
		while(q.size()) {
			int u = q.front();
			q.pop();
			for(int i = 0; i < 26; i++) {
				if(tr[u].next[i] != 0) {
					tr[tr[u].next[i]].fail = tr[tr[u].fail].next[i];
					q.push(tr[u].next[i]);
				} else {
					tr[u].next[i] = tr[tr[u].fail].next[i];
				}
			}
		} 
	}
	int work(string s) {
		getFail();
		int len = s.size();
		s = ' ' + s;
		int u = 0, ans = 0;
		for(int i = 1; i <= len; i++) {
			int x = s[i] - 'a';
			u = tr[u].next[x];
			int v = u;
			while(v && tr[v].cnt != -1) {
				ans += tr[v].cnt;
				tr[v].cnt = -1;
				v = tr[v].fail;
			}
		}
		return ans;
	}
};

signed main() {
	ACMachine USST;
	int n;
	cin >> n;
	USST.init();
	for(int i = 1; i <= n; i++) {
		string s;
		cin >> s;
		USST.build(s);
	}
	string s;
	cin >> s;
	cout << USST.work(s) << '\n';
}
```

## AC自动机 加强

![image-20231118024402942](C:\Users\10042\AppData\Roaming\Typora\typora-user-images\image-20231118024402942.png)

```cpp
#include <bits/stdc++.h>

using namespace std;
map<int, string> fmp;
map<int, int> Cnt;
struct ACMachine {
	int size;
	struct Node {
		int fail, cnt;
		int next[26]; 
		Node() {
			for(int i = 0; i < 26; i++) {
				next[i] = 0;
			}
			fail = 0;
			cnt = 0;
		}
	};
	vector<Node> tr;
	queue<int> q;
	void init() {
		fmp.clear();
		Cnt.clear();
		tr.clear();
		tr.push_back(Node());
		size = 0;
	}
	void build(string s) {
		int len = s.size(), u = 0;
		string t = s;
		s = ' ' + s;
		for(int i = 1; i <= len; i++) {
			int to = s[i] - 'a';
			if(!tr[u].next[to]) {
				tr.push_back(Node());
				size++;
				tr[u].next[to] = size;
			}
			u = tr[u].next[to];
		}
		fmp[u] = t; 
		tr[u].cnt++;
	}
	void getFail() {
		tr[0].fail = 0;
		for(int i = 0; i < 26; i++) {
			if(tr[0].next[i] != 0) {
				tr[tr[0].next[i]].fail = 0;
				q.push(tr[0].next[i]);
			}
		}
		while(q.size()) {
			int u = q.front();
			q.pop();
			for(int i = 0; i < 26; i++) {
				if(tr[u].next[i] != 0) {
					tr[tr[u].next[i]].fail = tr[tr[u].fail].next[i];
					q.push(tr[u].next[i]);
				} else {
					tr[u].next[i] = tr[tr[u].fail].next[i];
				}
			}
		} 
	}
	void work(string s) {
		getFail();
		int len = s.size();
		s = ' ' + s;
		int u = 0, ans = 0;
		for(int i = 1; i <= len; i++) {
			int x = s[i] - 'a';
			u = tr[u].next[x];
			int v = u;
			while(v) {
				if(tr[v].cnt) {
					Cnt[v]++;
				}
				v = tr[v].fail;
			}
		}
	}
};

void solve(int n) {
	ACMachine USST;
	USST.init();
	for(int i = 1; i <= n; i++) {
		string s;
		cin >> s;
		USST.build(s);
	}
	string s;
	cin >> s;
	USST.work(s);
	int ma = 0;
	for(auto [x, y] : Cnt) {
		ma = max(ma, y);
	}
	cout << ma << '\n';
	for(auto [x, y] : Cnt) {
		if(ma == y) {
			cout << fmp[x] << '\n';
		}
	}
}

signed main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	cout.tie(nullptr);
	int n;
	while(cin >> n, n) {
		solve(n);
	}
}
```

## AC自动机 二次加强

![image-20231118024500043](C:\Users\10042\AppData\Roaming\Typora\typora-user-images\image-20231118024500043.png)

```cpp
#include <bits/stdc++.h>

using namespace std;
map<int, string> to;
map<string, int> idx;
map<int, int> Cnt;
struct ACMachine {
	int size;
	struct Node {
		int fail, cnt;
		int next[26]; 
		Node() {
			for(int i = 0; i < 26; i++) {
				next[i] = 0;
			}
			fail = 0;
			cnt = 0;
		}
	};
	vector<Node> tr;
	queue<int> q;
	vector<int> que;
	void init() {
		Cnt.clear();
		tr.clear();
		tr.push_back(Node());
		que.push_back(0);
		size = 0;
	}
	void build(string s) {
		int len = s.size(), u = 0;
		string t = s;
		s = ' ' + s;
		for(int i = 1; i <= len; i++) {
			int to = s[i] - 'a';
			if(!tr[u].next[to]) {
				tr.push_back(Node());
				size++;
				tr[u].next[to] = size;
			}
			u = tr[u].next[to];
		}
		tr[u].cnt++;
		idx[t] = u;
	}
	void getFail() {
		tr[0].fail = 0;
		for(int i = 0; i < 26; i++) {
			if(tr[0].next[i] != 0) {
				tr[tr[0].next[i]].fail = 0;
				q.push(tr[0].next[i]);
				que.push_back(tr[0].next[i]);
			}
		}
		while(q.size()) {
			int u = q.front();
			q.pop();
			for(int i = 0; i < 26; i++) {
				if(tr[u].next[i] != 0) {
					tr[tr[u].next[i]].fail = tr[tr[u].fail].next[i];
					q.push(tr[u].next[i]);
					que.push_back(tr[u].next[i]);
				} else {
					tr[u].next[i] = tr[tr[u].fail].next[i];
				}
			}
		} 
	}
	void work(string s) {
		getFail();
		int len = s.size();
		s = ' ' + s;
		int u = 0, ans = 0;
		vector<int> du(size + 1);
		for(int i = 1; i <= len; i++) {
			int x = s[i] - 'a';
			u = tr[u].next[x];
			du[u]++;
		}
		// cout << size << ' ' << que.size() << '\n';
		for(int i = que.size() - 1; i >= 1; i--) {
			Cnt[que[i]] = du[que[i]];
			du[tr[que[i]].fail] += du[que[i]];
		}
	}
};

void solve(int n) {
	ACMachine USST;
	USST.init();
	for(int i = 1; i <= n; i++) {
		string s;
		cin >> s;
		USST.build(s);
		to[i] = s;
	}
	string s;
	cin >> s;
	USST.work(s);
	for(int i = 1; i <= n; i++) {
		cout << Cnt[idx[to[i]]] << '\n';
	}
}

signed main() {
	ios::sync_with_stdio(false);
	cin.tie(nullptr);
	cout.tie(nullptr);
	// freopen("1.in", "r", stdin);
	int n;
	cin >> n;		
	solve(n);
}
```

## sam

```CPP
struct SuffixAutomaton {
    static constexpr int ALPHABET_SIZE = 26, N = 1e5;
    struct Node {
        int len;
        int link;
        int next[ALPHABET_SIZE];
        Node() : len(0), link(0), next{} {}
    } t[2 * N];
    int cntNodes;
    SuffixAutomaton() {
        cntNodes = 1;
        std::fill(t[0].next, t[0].next + ALPHABET_SIZE, 1);
        t[0].len = -1;
    }
    int extend(int p, int c) {
        if (t[p].next[c]) {
            int q = t[p].next[c];
            if (t[q].len == t[p].len + 1)
                return q;
            int r = ++cntNodes;
            t[r].len = t[p].len + 1;
            t[r].link = t[q].link;
            std::copy(t[q].next, t[q].next + ALPHABET_SIZE, t[r].next);
            t[q].link = r;
            while (t[p].next[c] == q) {
                t[p].next[c] = r;
                p = t[p].link;
            }
            return r;
        }
        int cur = ++cntNodes;
        t[cur].len = t[p].len + 1;
        while (!t[p].next[c]) {
            t[p].next[c] = cur;
            p = t[p].link;
        }
        t[cur].link = extend(p, c);
        return cur;
    }
};
```

## PAM

```CPP
#include <bits/stdc++.h>
#define rep(i, a, b) for (int i=a; i<=b; i++)
#define drep(i, a, b) for (int i=a; i>=b; i--)
#define inf 1e9
using namespace std;
typedef long long ll;
const int maxn=300010;
char s[maxn];
int n;
struct PrefixAutomaton {
    int last;
    struct Node {
        int cnt, lenn, fail, son[27];
        Node(int lenn, int fail):lenn(lenn), fail(fail), cnt(0){
            memset(son, 0, sizeof(son));
        };
    };
    vector<Node> st;
    inline int newnode(int lenn, int fail=0) {
        st.emplace_back(lenn, fail);
        return st.size()-1;
    }
    inline int getfail(int x, int n) {
        while (s[n-st[x].lenn-1] != s[n]) x=st[x].fail;
        return x;
    }
    inline void extend(int c, int i) {
        int cur=getfail(last, i);
        if (!st[cur].son[c]) {
            int nw=newnode(st[cur].lenn+2, st[getfail(st[cur].fail, i)].son[c]);
            st[cur].son[c]=nw;
        }
        st[ last=st[cur].son[c] ].cnt++;
    }
    void init() {
        scanf("%s", s+1);
        n=strlen(s+1);
        s[0]=0;
        newnode(0, 1), newnode(-1);
        last=0;
        rep(i, 1, n) extend(s[i]-'a', i);
    }
    ll count() {
        drep(i, st.size()-1, 0) st[st[i].fail].cnt+=st[i].cnt;
        ll ans=0;
        rep(i, 2, st.size()-1) ans=max(ans, 1LL*st[i].lenn*st[i].cnt);
        return ans;
    }
}T;
int main() {
    T.init();
    printf("%lld\n", T.count());
}
```

## manacher

```cpp
std::vector<int> manacher(std::string s) {
    std::string t = "#";
    for (auto c : s) {
        t += c;
        t += '#';
    }
    int n = t.size();
    std::vector<int> r(n);
    for (int i = 0, j = 0; i < n; i++) {
        if (2 * j - i >= 0 && j + r[j] > i) {
            r[i] = std::min(r[2 * j - i], j + r[j] - i);
        }
        while (i - r[i] >= 0 && i + r[i] < n && t[i - r[i]] == t[i + r[i]]) {
            r[i] += 1;
        }
        if (i + r[i] > j + r[j]) {
            j = i;
        }
    }
    return r;
}
```

## kmp

```cpp
function<vector<int>(string)> kmp = [&](string s) {
     int n = s.size();
     s = ' ' + s;
     vector<int> nxt(n + 1);
     int j = 0;
     for(int i = 2, j = 0; i <= n; i++ ) {
          while(j && s[j + 1] != s[i]) {
               j = nxt[j];
          }
          if(s[j + 1] == s[i]) {
               j++;
          }
          nxt[i] = j;
     }
     return nxt;
};
```

## SA

```cpp
struct SA {
	string s;
	int n, m, len;
	vector<int> sa, rk, cnt;
	vector<int> height;
	vector<vector<int> > st;
	SA(string s) {
		n = s.size();
		s = ' ' + s;
		this -> s = s;
		len = 1;
		m = 128;
		while(len <= n) {
			len *= 2;
		}
		sa.assign(len + 1, 0);
		rk.assign(len + 1, 0);
		height.assign(len + 1, 0);
	}
	void getSA() {
		for(int i = 1; i <= n; i++) {
			sa[i] = i;
			rk[i] = s[i];
		}
		sort(sa.begin() + 1, sa.begin() + 1 + n, [&](int x, int y) {
			return rk[x] < rk[y];
		});
		int p = 0;
		for(int w = 1; w < n; w *= 2, m = p) {
			cnt.assign(m + 1, 0);
			auto id = sa;
			p = 0;
			for(int i = n; i > n - w; i--) {
				id[++p] = i;
			}
			for(int i = 1; i <= n; i++) {
				if(sa[i] > w) {
					id[++p] = sa[i] - w;
				}
			}
			for(int i = 1; i <= n; i++) {
				cnt[rk[id[i]]]++;
			}
			for(int i = 1; i <= m; i++) {
				cnt[i] += cnt[i - 1];
			}
			for(int i = n; i >= 1; i--) {
				sa[cnt[rk[id[i]]] --] = id[i];
			}
			auto copy = rk;
			int cnt = 0;
			for(int i = 1; i <= n; i++) {
				if(copy[sa[i]] == copy[sa[i - 1]] && copy[sa[i] + w] == copy[sa[i - 1] + w]) {
					rk[sa[i]] = cnt;
				} else {
					rk[sa[i]] = ++cnt;
				}
			}
		}
	}
	void getHeight() {
		for (int i = 1, k = 0; i <= n; ++i) {
		  	if (rk[i] == 0) continue;
		  	if (k) --k;
		  	while (s[i + k] == s[sa[rk[i] - 1] + k]) ++k;
		  	height[rk[i]] = k;
		}
	}
    void getST() {
		st.assign(n + 1, vector<int> (40));
		for(int i = 1; i <= n; i++) {
			st[i][0] = height[i];
		}
		for(int i = 1; i <= 30; i++) {
			int step = 1ll << (i - 1);
			for(int j = 1; j + step <= n; j++) {
				st[j][i] = min(st[j][i - 1], st[j + step][i - 1]);
			}
		}
	}
    /*
	int lcp(int x, int y) {
		if(y < x) {
			swap(x, y);
		}
		x++;
		int step = y - x + 1, mi = 1e18;
		for(int i = 30; i >= 0; i--) {
			int j = (1ll << i);
			if(step >= j) {
				step -= j;
				mi = min(mi, st[x][i]);
				x += j;
			}
		}
		return mi;
	} 倍增
	*/ 
    int lcp(int x, int y) { // O(1)
		if(y < x) {
			swap(x, y);
		}
		x++;
		int j = log2(y - x + 1);
		return min(st[x][j], st[y - (1ll << j) + 1][j]);
	}
	int LCP(int x, int y) {
		if(y < x) {
			swap(x, y);
		}
		return lcp(rk[x], rk[y]);
	}
	void work() {
		getSA();
		getHeight();
		getST();
	}
};
```



## mincut

最小割之算法模板
AcWing 2173. Dinic/ISAP求最小割
最小割之直接应用
AcWing 2279. 网络战争
01分数规划，对答案进行二分查找，不断建图并使用最小割来作二分查找的check函数
最小化的是：选定边集的总权值和 - 边集的大小*枚举的平均值，则合并一下，最小化（选定的边集中所有边权减平均值的权值和）即可，注意减成负数的边必选，直接累计其原始边权值（不是减去平均值后的负数），并将这正反一对边的容量置零，防止跑最小割时出错

AcWing 2280. 最优标号
将多个数按二进制位处理，2^31-1，即31位，根据当前数在枚举的二进制位上，是0或1来建图，建31次图跑31次最小割

插曲：最小路径覆盖
最小路径覆盖：选择最少的路径（从任意一个点开始，到任意一个不为起点的点结束），将所有点都覆盖
最小路径覆盖 = 总点数 - 最大匹配

P2764 最小路径覆盖问题
将所有点都拆成一对出点和入点，源点连所有入点，所有出点连汇点，容量均为1
对于给出的图中边（a，b），我们残留网络连边（a的入点，b的出点，INF），这样就构造了一条S->a(入)->b(出)->T的路径，跑最大流即可获的最大匹配数，即在每个点的出点入点都最多选择一次的情况下，选择了最大数量（即流量）的边，而对于每一条边，都可以将这一对点合并，则选择了多少条边就最多可以合并多少个点，每合并一个点就少一条路径，所以最小路径覆盖 = 总点数 - 最大匹配
最小割之最大权闭合图
AcWing 961. 最大获利

```cpp
int main()
{
    scanf("%d%d", &n, &m);

    for(int i = 1; i <= n; i++) scanf("%d", &p[i]), p[i] *= -1; //点权是要减去的，是负的

    memset(h, -1, sizeof h); //初始化邻接表

    S = 0, T = n + 1; //源点、汇点
    while(m--)
    {
        int a, b, c;
        scanf("%d%d%d", &a, &b, &c);

        deg[a] += c, deg[b] += c; //记录从每个点出发的边的权值之和
        add(a, b, c, c); //无向边（四条边合并成两条边）
    }
    
    /*
    void add(int a, int b, int c, int d) //添加边
    {
        e[idx] = b, w[idx] = c, ne[idx] = h[a], h[a] = idx++;
        e[idx] = a, w[idx] = d, ne[idx] = h[b], h[b] = idx++;
    }
    */

    //保证 U + 2 * g - 2 * p[i] - deg[i] >= 0，本题 g = 0，U >= 2 * p[i] + deg[i]
    int U = 0;
    for(int i = 1; i <= n; i++) U = max(U, 2 * p[i] + deg[i]);

    for(int i = 1; i <= n; i++)
    {
        add(S, i, U, 0); //从源点向每个点连一条容量是 U 的边
        add(i, T, U - 2 * p[i] - deg[i], 0); //从每个点向汇点连一条容量是 U - 2 * g - 2 * p[i] - deg[i] 的边
    }

    printf("%d\n", (U * n - dinic()) / 2);

    return 0;
}

```



最大权闭合图：选择若干个点构成点集，不存在从该点集内部指向点集外部的有向边（点集闭合），同时使选择的点总权值和最大

适用模型：对于一个点集A，选择其中一点可以得到相应点权的价值，但是对于属于点集A的任何一点a，可以选择它的前提是选择若干个它的必要点，并会损失选择这些必要点的点权的价值（若能选a的前提是选择b和c，则我们最终选择a能得到的价值为：-b-c+a），要求使最终价值最大化，和对应的一种方案。

建图（残留网络）：
我们称能获得价值的点为正权点，而选择它需要先选择的若干必要点位称为负权点
1、从源点向正权点连边，容量为正权值
2、从负权点向汇点连边，容量为权值绝对值
3、其他点与点之间，按给出的关系连边，容量为INF ==> 最小割一定不会选容量为INF的边，则保证了选出的点集一定是合理的（点集闭合）

答案：所有正权点权值和，减去最小割，

1、理解：
因不能使容量为INF的边成为割边，我们假设流量从正权点g流出，一定要经过所有它能到达的负权点，若不能把所有负权点流满或者刚刚好流满，则说明选这个正权点g能得到的代价，小于等于损失的代价，这个点就是废点，我们肯定不选，注意：此时边S->g的容量为0，若全流满，说明选这个正权点g能得到的代价，大于损失的代价，注意：此时边S->g的容量有剩余。

2、证明：
证明可得（此处省略），任何一个闭合图都可以对应一个简单割，而构造这个简单割就对应了我们不选哪些"废点"，和选了哪些"好点"对应的损耗。
除了中间因相互制约而出现的容量为INF的边以外，剩下的边有两种，从源点出来向正权点的边，和负权点通往汇点的边，且割边必须是这两种边：
①若割边是从源点出来的边，代表在最大流中，这条边比它的必要点的权值总和小，是废点，割边容量为该正权点的点权
②若割边是通往汇点的边，代表在最大流中，这条边比它的必要点的权值总和大，是好的，割边容量为该负权点的点权
综上，最小割容量表示的是：所有没有选上的正权点的权值，和所有负权点中，与选上了的正权点相关的点的权值，则正权点权值总和与最小割相减，自然就减去了没选的正权点，还减去了选了的正权点对应的负权消耗。

最大权闭合图的方案：注意，根据以上"1、理解"推导可知，最大权闭合子图选上的点是割集S，即从源点沿容量不为零的边能遍历到的所有点

最小割之最大密度子图
AcWing 2324. 生活的艰辛
最大化子图的密度，其密度的分式：边权和 / 点集大小
符合01分数规划，则使用二分来枚举答案，不断建新图和跑最小割来作二分查找的check（）函数
二分最大化的是E-g×V，则改为最小化g×V-E，就可以使用最小割，式子最终化简为 U×n - 最小割
（这里的g指的是二分枚举的密度值，n为点的总数）

建图（残留网络）：
1、根据题意点与点直接的关系建图，双向关系建容量相等的一对边
2、为保证边权始终为正，源点流出和流入汇点的边都加上一个偏移量U，此操作在这里的合理性可证，其他情况不可随意对图中的边权做加减操作
3、则源点向所有点连边，容量为U，所有点向汇点连边，容量为U + 2g - du【i】，这里的U是权值偏移量，g指的是二分枚举的密度值，du【i】指的是该点的度数（出度+入度）

最小割之最小点权覆盖集
AcWing 2325. 有向图破坏
最小点权覆盖集：指选择任意数量的点构成点集，并使图中任意一条边都至少有一个顶点被选中，最小化选中点的权值和

本题可转化为，求：对于图中任意一条边，都至少选中其中一个顶点，并最小化选中点的点权和
将所有点分为“出”和“入”两种，
源点向所有出点连边，容量为出点权，
所有入点向源点连边，容量为入点权，
然后对于原图中的有向边（X，Y），从X的出点向Y的入点连边，容量为INF。

①与最大权闭合图原理类似，最小割一定不会选容量为INF的边，即割边一定是S->X，或X->T这样的边，
②保证了对于一条容量INF的边，两个点不会同时被选，因为此时即存在了S->X->Y->T这样的边，不符合割的定义，
③对于割边，若它是起点为S的边，则选择这个入点，若它是终点为T的边，选择这个出点。

最小点权覆盖集的值 等于 最小割的容量 的充分必要性证明：
1、一个简单割是否是一个点覆盖集？对于原图的边（X，Y），在残留网络中以S->X->Y->T这样的形式存在，X-Y的边必不是割边，若S->X，Y->T也均不是割边的话，即S->X->Y->T这三条边都在割的左部分或右部分，但S和T不应在同一部分。所以S->X和Y->T必有其一且只能是其一被选（②），即满足点覆盖集的概念。
2、一个点覆盖集是否是一个简单割？对于边（X，Y），若选择了X，则将S->X设为割边，若选择了Y，则将Y->T设为割边，证明类似，略。

最小点权覆盖集的方案：在残留网络中，若存在容量不为0，且是割集S->割集T的边，则其为割边，为S->Z的类型则选择出点Z，为Z->T的类型则选择入点Z。

最小割之最大点权独立集
AcWing 2326. 王者之剑
独立集：指在选出的点集中，任意两点间都不存在任何有向边
点覆盖集的补集一定是一个独立集，反之也成立：独立集的补集一定是一个点覆盖集

可证（证明略），任意一个独立集都可以对应一个覆盖集，那么公式：
最大点权独立集 = 所有点权总和 - 最小点权覆盖集

适用模型：对于图中的给定的二元关系，它们是互斥的，即选择了A就不能选择B，那么我们就连边A->B（容量为INF）表示这两个点应该相互独立，求选择若干个独立点的最大价值

建图（残留网络）：
创建虚拟源汇点S，T，将所有方格（x，y）按（x+y）%2，分成奇偶两类，一类全部由源点进来，一类全部到达汇点，容量均为点的权值

那么在问题转化为求最小点权覆盖集后，很明显，同一有向边的两个端点不能同时选，对应到问题中就是当前位置与上下左右四个方位相邻的位置连边，容量为INF，表示不能同时选择任何两个相邻的方格。
————————————————
