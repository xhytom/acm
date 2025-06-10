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
const ll N = 500005;
//head

struct Node{
    int x;
};

Node tr[4 * N];

void build(int u, int L, int R)
{
    tr[u].x = 0;
    if(L == R) return;
    build(2 * u, L, (L + R) / 2);
    build(2 * u + 1,(L + R) / 2 + 1, R);
}

void pushup(int u, int y)
{
    //tr[u/2].x += y;
}

void change(int u, int L, int R, int x, int y)
{
    tr[u].x += y;
    if(L == R)
    {
        return;
    }else{
        if(x <= (L + R) / 2)
        {
            change(2 * u, L, (L + R) / 2, x, y);
        }else{
            change(2 * u + 1, (L + R) / 2 + 1, R, x, y);
        }
    }
    pushup(u, y);
}

ll query(int u, int L, int R, int l, int r)
{
    if(l <= L && r >= R)
    {
        return tr[u].x;
    }
    if(r <= (L+R)/2 )
    {
        return query(2*u, L, (L+R)/2, l, r);
    }
    if(l >= (L+R)/2+1)
    {
        return query(2*u + 1, (L+R)/2+ 1, R, l, r);
    }
    return query(2*u, L, (L+R)/2, l, r) + query(2*u + 1, (L+R)/2+ 1, R, l, r);
}
int a[N];

void solve()
{
    int n;
    cin >> n;
    build(1, 1, n);
    for(int i = 1 ; i <= n ; i++ )
    {
        cin >> a[i];
    }
    for(int i = 1 ; i <= n ; i++ )
    {
        change(1, 1, n, i, a[i]);
    }
    while(1)
    {
        string s;
        cin >> s;
        if(s[0] == 'E')
        {
            break;
        }
        if(s[0] == 'A')
        {
            int x, y;
            cin >> x >> y;
            change(1, 1, n, x, y);
        }
        if(s[0] == 'S')
        {
            int x, y;
            cin >> x >> y;
            change(1, 1, n, x, -y);
        }
        if(s[0] == 'Q')
        {
            int x, y;
            cin >> x >> y;
            cout << query(1, 1, n, x, y) << "\n";
        }
    }
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    int t;
    cin >> t;
    //t = 1;
    for(int i = 1 ; i <= t ; i++ )
    {
        cout << "Case " << i << ":\n";
        solve();
    }
    return 0;
}