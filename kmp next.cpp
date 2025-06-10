#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
const ll MOD = 998244353;
const ll N = 200005;
//head

int nxt[N];
char s[N];

void solve()
{   
    cin >> s + 1;
    int n = strlen(s + 1);
    nxt[1] = 0;
    int j = 0;
    for(int i = 2 ; i <= n ; i++ )
    {
        while(j > 0 && s[j + 1] != s[i])
        {
            j = nxt[j];
        }
        if(s[j + 1] == s[i])
        {
            j++;
        }
        nxt[i] = j;
    }
    for(int i = 1 ; i <= n ; i++ )
    {
        cout << nxt[i] << " \n"[i == n];
    }

}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    int t;
    t = 1;
    //cin >> t;
    while(t--)
    {
        solve();
    }
    return 0;
}   