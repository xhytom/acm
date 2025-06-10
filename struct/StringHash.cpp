// Problem: 3/2 Square Strings
// Contest: NowCoder
// URL: https://ac.nowcoder.com/acm/contest/85132/G
// Memory Limit: 2097152 MB
// Time Limit: 8000 ms
// 
// Powered by CP Editor (https://cpeditor.org)

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
#pragma GCC optimize("O3,unroll-loops") 
//这行告诉GCC编译器使用O3优化级别和循环展开。O3是GCC提供的最高优化级别，它会尝试使用所有的程序优化策略。"unroll-loops"是一个特定的优化选项，它会尝试将循环展开以减少循环的开销。

#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt") 

//这行告诉GCC编译器生成的代码应该针对支持AVX2，BMI，BMI2，LZCNT和POPCNT指令集的CPU。这些都是特定的CPU指令集，可以提高代码的性能，但是生成的代码可能无法在不支持这些指令集的CPU上运行。
#pragma GCC optimize("Ofast")
#pragma GCC target("avx", "sse2")
#pragma GCC optimize("inline")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-fwhole-program")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-fstrict-overflow")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-skip-blocks")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("-funsafe-loop-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
using i64 = long long;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
#define debug(x) cerr << #x << " = " << (x) << endl;
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

using u64 = uint64_t;
struct LongestCommonPrefix {
    int n;
    vector<int> p, rank;
    vector<vector<int>> st;
    LongestCommonPrefix(const string &s) : n(s.size()), p(n), rank(n) {
        int k = 0;
        vector<int> q, count;

        for (int i = 0; i < n; i += 1)
            p[i] = i;

        sort(p.begin(), p.end(), [&](int i, int j) {
            return s[i] < s[j];
        });

        for (int i = 0; i < n; i += 1)
            rank[p[i]] = i and s[p[i]] == s[p[i - 1]] ? rank[p[i - 1]] : k++;

        for (int m = 1; m < n; m *= 2) {
            q.resize(m);

            for (int i = 0; i < m; i += 1)
                q[i] = n - m + i;

            for (int i : p)
                if (i >= m)
                    q.push_back(i - m);

            count.assign(k, 0);

            for (int i : rank)
                count[i] += 1;

            for (int i = 1; i < k; i += 1)
                count[i] += count[i - 1];

            for (int i = n - 1; i >= 0; i -= 1)
                p[count[rank[q[i]]] -= 1] = q[i];

            auto cur = rank;
            cur.resize(2 * n, -1);
            k = 0;

            for (int i = 0; i < n; i += 1)
                rank[p[i]] = i and cur[p[i]] == cur[p[i - 1]] and
                             cur[p[i] + m] == cur[p[i - 1] + m]
                             ? rank[p[i - 1]]
                             : k++;
        }

        st.emplace_back(n);

        for (int i = 0, k = 0; i < n; i += 1) {
            if (not rank[i])
                continue;

            k = max(k - 1, 0);
            int j = p[rank[i] - 1];

            while (i + k < n and j + k < n and s[i + k] == s[j + k])
                k += 1;

            st[0][rank[i]] = k;
        }

        for (int i = 1; (1 << i) < n; i += 1) {
            st.emplace_back(n - (1 << i) + 1);

            for (int j = 0; j <= n - (1 << i); j += 1)
                st[i][j] = min(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
        }
    }
    int get(int i, int j) {
        if (i == j)
            return n - i;

        if (i == n or j == n)
            return 0;

        i = rank[i];
        j = rank[j];

        if (i > j)
            swap(i, j);

        int k = 64 - __builtin_clzll(u64(j - i)) - 1;
        return min(st[k][i + 1], st[k][j - (1 << k) + 1]);
    }
    
    int rankget(int i, int j) {
        if (i == j)
            return 1000000000LL;

       
        if (i > j)
            swap(i, j);

        int k = 64 - __builtin_clzll(u64(j - i)) - 1;
        return min(st[k][i + 1], st[k][j - (1 << k) + 1]);
    }
};

vector<tuple<int, int, int>> run(const string &s) {
    int n = s.size();
    auto r = s;
    reverse(r.begin(), r.end());
    LongestCommonPrefix lcp(s), lcs(r);
    vector<tuple<int, int, int>> runs;

    for (bool inv : {
                false, true
            }) {
        vector<int> lyn(n, n), stack;

        for (int i = 0; i < n; i += 1) {
            while (not stack.empty()) {
                int j = stack.back(), k = lcp.get(i, j);

                if (i + k < n and ((s[i + k] > s[j + k]) ^ inv))
                    break;

                lyn[j] = i;
                stack.pop_back();
            }

            stack.push_back(i);
        }

        for (int i = 0; i < n; i += 1) {
            int j = lyn[i], t = j - i, l = i - lcs.get(n - i, n - j),
                r = j + lcp.get(i, j);

            if (r - l >= 2 * t)
                runs.emplace_back(t, l, r);
        }
    }
    sort(runs.begin(), runs.end());
    runs.erase(unique(runs.begin(), runs.end()), runs.end());
    return runs;
}




struct HashValue {
    int x = 0, y = 0;
    HashValue(int _x) {
        x = y = _x;
    } 
    HashValue(int _x, int _y) {
        x = _x;
        y = _y;
    }
    HashValue(){}
    
    friend HashValue operator+(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x + rhs.x, lhs.y + rhs.y};
    }

    friend HashValue operator%(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x % rhs.x, lhs.y % rhs.y};
    }

    friend HashValue operator*(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x * rhs.x, lhs.y * rhs.y};
    }

    friend HashValue operator-(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x - rhs.x, lhs.y - rhs.y};
    }

    friend std::ostream &operator<<(std::ostream &os, HashValue value) {
        return os << value.x << " " << value.y;
    }

    friend bool operator==(const HashValue &lhs, const HashValue &rhs) {
        return lhs.x == rhs.x && lhs.y == rhs.y;
    }

    friend bool operator<(const HashValue &lhs, const HashValue &rhs) {
        return lhs.x == rhs.x ? lhs.y < rhs.y : lhs.x < rhs.x;
    }
};



const struct HashBase {
    std::vector<HashValue> fact;
    HashValue P, base;
    int nextPrime(int n) {
        int flag = 1;
        while (flag) {
            n++;
            flag = 0;
            for (int i = 2; i * i <= n; i++) {
                if (n % i == 0) {
                    flag = 1;
                    break;
                }
            }
        }
        return n;
    }
    HashBase(int n) : fact(n + 5){
        const int P1 = nextPrime(rnd(900000000, 1000000000));
        const int P2 = nextPrime(rnd(900000000, 1000000000));
        P = HashValue(P1, P2);
        const int base1 = nextPrime(rnd(200, 10000));
        const int base2 = nextPrime(rnd(200, 10000));
        base = HashValue(base1, base2);
        fact[0] = HashValue(1);
        for (int i = 1; i <= n + 3; i++) {
            fact[i] = fact[i - 1] * base % P;
        }
    }
    HashBase(){}    
}B(N);

struct DynamicStringHash {
    std::vector<HashValue> pre, bac;
    DynamicStringHash() : pre(1){
        pre.reserve(N);
        pre[0] = {0, 0};
    }
    int size() {
        return pre.size();
    }
    void push_back(char ch) {
        pre.push_back(pre.back() * B.base % B.P + HashValue(ch) % B.P) ;
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return ((pre[r] - (pre[l - 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + get(l, r)) % B.P;
    }
};

struct StringHash {
    std::vector<HashValue> pre, bac;
    StringHash(std::string s) : pre(s.size() + 1), bac(s.size() + 1){
        pre[0] = {0, 0};
        for (int i = 1; i < s.size(); i++) {
            pre[i] = (pre[i - 1] * B.base % B.P + HashValue(s[i])) % B.P;
        }
        bac[s.size()] = {0, 0};
        for (int i = s.size() - 1; i >= 0; i--) {
            bac[i] = (bac[i + 1] * B.base % B.P + HashValue(s[i])) % B.P;
        }

    }
    HashValue getrev(int l, int r) {
        if (l > r) return HashValue();
        return ((bac[l] - (bac[r + 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return ((pre[r] - (pre[l - 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + get(l, r)) % B.P;
    }
    HashValue concatrev(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + getrev(l, r)) % B.P;
    }
};

signed main()
{  
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    std::cout << std::fixed << std::setprecision(10);
    std::string s, t;
    std::cin >> s >> t;
    int sn = s.size(), tn = t.size();
    auto g = '#' + s + '$' + t;
    t = '#' + t;
    auto runs = run(t);
    StringHash T(t);
    std::map<HashValue, int> cnt;
    std::map<HashValue, std::pair<int, int>> pos;
    for (auto [P, l, r] : runs) {
        r--;
        for (int p = P; p <= (r - l + 1); p += P) {
            for (int i = l; i < l + p; i++) {
                if (i + 2 * p - 1 > r) break;
                if (i + 2 * p - 1 <= r) {
                    int x = (r - i + 1) / p - 1;
                    HashValue val = T.get(i, i + p - 1);
                    cnt[val] += x;
                    pos[val] = {i, i + p - 1};
                }
            }
        }
    }
    
    // 1 ~ sn  sn + 1,   sn + 2 ~ sn + tn + 1
    LongestCommonPrefix sa(g);
    i64 ans = 0;
    
    std::vector<int> pre(g.size() + 2);
    for (int i = 1; i < g.size(); i++) {
        if (1 <= sa.p[i] && sa.p[i] <= sn) {
            pre[i] = 1;
        }
        pre[i] += pre[i - 1];
        // for (int j = sa.p[i]; j < g.size(); j++) {
            // std::cerr << g[j];
        // }
        // std::cerr << "\n";
        // std::cerr << sa.p[i] << " ";
    }
    // rank[i] 表示 第 i 个后缀的排名
    // p[i] 表示第 i 名的字符串的位置
    
    
    for (auto [val, count] : cnt) {
        auto [l, r] = pos[val];
        int len = r - l + 1;
        int pos = sa.rank[l + sn + 1];
        int lo = 1, hi = pos;
        while (lo < hi) {
            int mid = lo + hi >> 1;
            if (sa.rankget(mid, pos) >= len) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        int LEFT = lo;
        lo = pos, hi = g.size() - 1;
        while (lo < hi) {
            int mid = lo + hi + 1 >> 1;
            if (sa.rankget(pos, mid) >= len) {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        int RIGHT = hi;
        // test(pos, l, r, LEFT, RIGHT, count);
        ans += (1LL * pre[RIGHT] - pre[LEFT - 1]) * count % MOD;
        ans %= MOD;
    }
    std::cout << ans << "\n";
    
    
    return 0;
}