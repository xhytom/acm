#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
const ll MOD = 998244353;
//const ll MOD = 1e9+7;
const ll N = 200005;
//head

struct Node{
    ll sum, lazy, size;
};
Node tr[N * 4];
ll a[N];

void pushup(int u, int L, int R)
{
    tr[u].sum = tr[u << 1].sum + tr[u << 1 | 1].sum;
}


void build(int u, int L, int R)
{
    int mid = L + R >> 1;
    tr[u].size = R - L + 1;
    tr[u].sum = tr[u].lazy = 0;
    if(L == R)
    {
        tr[u].sum = a[L];
        return;
    }
    build(u << 1, L, mid);
    build(u << 1 | 1, mid + 1, R);
    pushup(u, L, R);

}

void pushdown(int u)
{
    auto &root = tr[u], &left = tr[u << 1], &right = tr[u << 1 | 1]; 
    if(root.lazy)
    {
        left.sum += root.lazy * left.size;
        left.lazy += root.lazy;
        right.sum += root.lazy * right.size;
        right.lazy += root.lazy;
        root.lazy = 0;
    }
}

void pushup(int u)
{
    tr[u].sum = tr[u << 1].sum + tr[u << 1 | 1].sum;
}

ll query(int u, int L, int R, int l, int r)
{
    int mid = L + R >> 1;
    if(l <= L && R <= r)
    {
        return tr[u].sum;
    }
    ll ans = 0;
    pushdown(u);
    if(l <= mid)
    {
        ans += query(u << 1, L, mid, l, r);
    }
    if(r > mid)
    {
        ans += query(u << 1 | 1, mid + 1, R, l, r);
    }
    return ans;
}

void modify(int u, int L, int R, int l, int r, int x)
{
    int mid = L + R >> 1;
    if(l <= L && R <= r)
    {
        tr[u].lazy += x;
        tr[u].sum += x * tr[u].size;
        return;
    }
    pushdown(u);
    if(l <= mid)
    {
        modify(u << 1, L, mid, l , r, x);
    }
    if(r > mid)
    {
        modify(u << 1 | 1, mid + 1, R, l, r, x);
    }
    pushup(u);
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    int n, m;
    cin >> n >> m;
    for(int i = 1 ; i <= n ; i++ )
    {
        cin >> a[i];
    }
    build(1, 1, n);
    for(int i = 1 ; i <= m ; i++ )
    {
        int op;
        cin >> op;
        if(op == 1)
        {
            int x, y, d;
            cin >> x >> y >> d;
            modify(1, 1, n, x, y, d);
        }else{
            int x, y;
            cin >> x >> y;
            cout << query(1, 1, n, x, y) << "\n";
        }
        for(int i = 1 ; i <= n ; i++ )
            cout << query(1, 1, n, i, i) << " \n"[i == n];
    }
    return 0;
}   