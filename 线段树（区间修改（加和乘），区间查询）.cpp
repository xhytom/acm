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
ll P = 571373;
//const ll P = 1e9+7;
const ll N = 100005;
//head

struct Node{
    ll sum, mul, add, size;
} tr[4 * N];
ll a[N];

void pushup(int u)
{
    tr[u].sum = (tr[u << 1].sum % P + tr[u << 1 | 1].sum % P) % P;
}

void pushdown(int u)
{
    auto &root = tr[u], &left = tr[u << 1], &right = tr[u << 1 | 1];
    root.mul  %= P, root.add %= P;
    left.sum  *= root.mul;               left.sum  %= P;
    left.sum  += root.add * left.size;   left.sum  %= P;
    right.sum *= root.mul;               right.sum %= P;
    right.sum += root.add * right.size;  right.sum %= P;
    left.add  *= root.mul;               left.add  %= P;
    left.mul  *= root.mul;               left.mul  %= P;
    right.add *= root.mul;               right.add %= P;
    right.mul *= root.mul;               right.mul %= P;
    left.add  += root.add;               left.add  %= P;
    right.add += root.add;               right.add %= P;
    root.mul  = 1;               
    root.add  = 0;
}

void build(int u, int L, int R)
{
    int mid = L + R >> 1;
    tr[u].size = R - L + 1;
    tr[u].mul = 1;
    tr[u].add = 0;
    if(L == R)
    {
        tr[u].sum = a[L] % P;
        return;
    }
    build(u << 1, L, mid);
    build(u << 1 | 1, mid + 1, R);
    pushup(u);
}

void modify_add(int u, int L, int R, int l, int r, int x)
{
    int mid = L + R >> 1;
    if(l <= L && R <= r)
    {
        tr[u].sum += tr[u].size * x;    tr[u].sum %= P;

        tr[u].add += x;                 tr[u].add %= P;
        return;
    }
    pushdown(u);
    if(l <= mid)
    {
        modify_add(u << 1, L, mid, l, r, x);
    }
    if(r >= mid + 1)
    {
        modify_add(u << 1 | 1, mid + 1, R, l, r, x);
    }
    pushup(u);
}

void modify_mul(int u, int L, int R, int l, int r, int x)
{
    int mid = L + R >> 1;
    if(l <= L && R <= r)
    {
        tr[u].sum *= x; tr[u].sum %= P;
        tr[u].add *= x; tr[u].add %= P;
        tr[u].mul *= x; tr[u].mul %= P;
        return;
    }
    pushdown(u);
    if(l <= mid)
    {
        modify_mul(u << 1, L, mid, l, r, x);
    }
    if(r >= mid + 1)
    {
        modify_mul(u << 1 | 1, mid + 1, R, l, r, x);
    }
    pushup(u);
}

ll query(int u, int L, int R, int l, int r)
{
    if(l <= L && R <= r)
    {
        return tr[u].sum % P;
    }
    pushdown(u);
    ll ans = 0;
    int mid = L + R >> 1;
    if(l <= mid)
    {
        ans += query(u << 1, L, mid, l, r);
        ans %= P;
    }
    if(r >= mid + 1)
    {
        ans += query(u << 1 | 1, mid + 1, R, l, r);
        ans %= P;
    }
    pushup(u);
    return ans % P;
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    //freopen("my.out", "w", stdout);
    int n, m;
    cin >> n >> m >> P;
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
            int x, y, k;
            cin >> x >> y >> k;
            k %= P; 
            modify_mul(1, 1, n, x, y, k);
        }else if (op == 2)
        {
            int x, y, k;
            cin >> x >> y >> k;
            k %= P;
            modify_add(1, 1, n, x, y, k);
        }else{
            int x, y;
            cin >> x >> y;
            cout << query(1, 1, n, x, y) % P<< "\n";
        }
        /*for(int i = 1 ; i <= n ; i++ )
        {
            cout << query(1, 1, n, i, i) << " \n"[i == n];
        }*/
    }
    return 0;
}