#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
const ll MOD = 998244353;
//const ll MOD = 1e9+7;
const ll N = 100005;
//head

inline int read()
{
    int x=0,f=1;char ch=getchar();
    while (ch<'0'||ch>'9'){if (ch=='-') f=-1;ch=getchar();}
    while (ch>='0'&&ch<='9'){x=x*10+ch-48;ch=getchar();}
    return x*f;
}

ll a[N];
ll f[25][N];//0 - 1   1 - 2  2 - 4

int main()
{  
    //fastio
    //freopen("1.in","r",stdin);
    int n = read();
    for(int i = 1 ; i <= n ; i++ )
    {
        a[i] = read();
        f[0][i] = a[i];
    }
    for(int i = 1 ; i <= 22 ; i++ )
    {
        for(int j = 1 ; j + (1 << i) - 1 <= n ; j++ )
        {
            f[i][j] = max(f[i-1][j], f[i-1][j + (1 << i - 1)]);
        }
    }
    auto query = [&](int l, int r) {
        int len = __lg(r - l + 1);
        return(min(f[len][l], f[len][r - (1 << len) + 1]));
    };

    return 0;
}