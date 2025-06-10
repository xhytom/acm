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
using i64 = long long;
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
struct Node {
    Node *l = nullptr;
    Node *r = nullptr;
    int val = 0;
    int sz = 1;
    int cnt = 1;
    int prio = 0;
    bool tag = 0;
    Node(){}
    Node (int _val) : val(_val) {
        prio = mrand();

    }
    void pushup() {
        sz = cnt;
        if (l != nullptr) sz += l -> sz;
        if (r != nullptr) sz += r -> sz;
    }

    void pushdown() {
        if (tag) {
            std::swap(l, r);
            l -> tag ^= 1;
            r -> tag ^= 1;
            tag = false;
        }
    }


} pool[21 << 10];

struct Treap{
    int idx = 0;

    Node* root = nullptr;
    Node* newnode(int x) {
        pool[++idx] = Node(x);
        return &pool[idx];
    }

    int siz(Node *cur) {
        return cur == nullptr ? 0 : cur -> sz;
    }


    std::pair<Node*, Node*> split(Node *cur, int sz) {
        if (cur == nullptr) return {nullptr, nullptr};
        cur -> pushdown();
        if (sz <= siz(cur -> l)) {
            auto tmp = split(cur -> l, sz);
            cur -> l = tmp.second;
            cur -> pushup();
            return {tmp.first, cur};
        } else {
            auto tmp = split(cur -> r, sz - siz(cur -> l) - 1);
            cur -> r = tmp.first;
            cur -> pushup();
            return {cur, tmp.second};
        }
    }

    Node *merge(Node *sm, Node *bg) {
        if (sm == nullptr && bg == nullptr) return nullptr;
        if (sm == nullptr && bg != nullptr) return bg;
        if (sm != nullptr && bg == nullptr) return sm;
        sm -> pushdown();
        bg -> pushdown();
        if (sm -> prio < bg -> prio) {
            sm -> r = merge(sm -> r, bg);
            sm -> pushup();
            return sm;
        } else {
            bg -> l = merge(sm, bg -> l);
            bg -> pushup();
            return bg;
        }
    }

    void insert(int x) {
        debug(root);
        auto tmp = split(root, x);
        auto left = split(tmp.first, x - 1);
        Node *cur;
        if (left.second == nullptr) cur = newnode(x);
        Node *p = merge(left.first, left.second == nullptr ? cur : left.second);
        root = merge(p, tmp.second);
    }

    void reverse(int l, int r) {
        auto p1 = split(root, l - 1);
        auto p2 = split(p1.second, r - l + 1);
        p2.first -> tag ^= 1;
        root = merge(p1.first, merge(p2.first, p2.second));
    }

    void print(Node* cur) {
        if (cur -> l != nullptr) {
            print(cur -> l);
        }
        std::cout << cur -> val << " ";
        if (cur -> r != nullptr) {
            print(cur -> r);
        }
    }

} T;

signed main()
{  
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    int n;
    std::cin >> n;
    for (int i = 1; i <= n; i++) {
        T.insert(i);
        debug(i);
    }
    int q;
    std::cin >> q;
    while (q--) {
        int l, r;
        std::cin >> l >> r;
        T.reverse(l, r);
    }
    T.print(T.root);



    return 0;
}