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
    int minx, cntminx;
};

ll a[N];

Node tr[4 * N];

void pushup(int u, int L, int R)
{
    if(tr[u << 1].minx < tr[u << 1 | 1].minx)
    {
        tr[u].minx = tr[u << 1].minx;
        tr[u].cntminx = tr[u << 1].cntminx;
    }
    if(tr[u << 1].minx > tr[u << 1 | 1].minx)
    {
        tr[u].minx = tr[u << 1 | 1].minx;
        tr[u].cntminx = tr[u << 1 | 1].cntminx;
    }
    if(tr[u << 1].minx == tr[u << 1 | 1].minx)
    {
        tr[u].minx = tr[u << 1 | 1].minx;
        tr[u].cntminx = tr[u << 1].cntminx + tr[u << 1 | 1].cntminx;
    }
}


void build(int u, int L, int R)
{
    int mid = L + R >> 1;
    if(L == R)
    {
        tr[u].minx = a[L];
        tr[u].cntminx = 1;
        return;
    }
    build(u << 1, L, mid);
    build(u << 1 | 1, mid + 1, R);
    pushup(u, L, R);

}

void change(int u, int L, int R, int x, int y)
{
    int mid = L + R >> 1;
    if(L == R)
    {
        tr[u].minx += y;
        return;
    }
    if(x <= mid)
    {
        change(u << 1, L, mid, x, y);
    }
    if(x > mid)
    {
        change(u << 1 | 1, mid + 1, R, x, y);
    }
    pushup(u, L, R);
}

pair<int, int> query(int u, int L, int R, int l, int r)
{
    int mid = L + R >> 1;
    if(l <= L && R <= r)
    {
        return {tr[u].minx, tr[u].cntminx};
    }
    if(r <= mid)
    {
        return query(u << 1, L, mid, l, r);
    }
    if(l >= mid + 1)
    {
        return query(u << 1 | 1, mid + 1, R, l, r);
    }
    auto s1 = query(u << 1, L, mid, l, r);
    auto s2 = query(u << 1 | 1, mid + 1, R, l, r);
    if(s1.first < s2.first)
    {
        return s1;
    }
    if(s1.first > s2.first)
    {
        return s2;
    }
    return {s1.first, s1.second + s2.second};
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
        int op, x, y;
        cin >> op >> x >> y;
        if(op == 1)
        {
            change(1, 1, n, x, y);
        }else{
            auto [_,__] =  query(1, 1, n, x, y);
            cout << _ << " " << __ << "\n";
        }
    }
    return 0;
}