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

int MAXN = 50;
int p[N], pr[N], idx;

void build()
{
    for(int i = 2 ; i < MAXN ; i++ )
    {
        if(!p[i])
        {
            p[i] = i;
            pr[++idx] = i;
        }
        for(int j = 1 ; j <= idx && pr[j] * i < MAXN ; j++ )
        {
            p[i * pr[j]] = pr[j];
            if(p[i] == pr[j]) break;
        }
    }
    /*for(int i = 1 ; i <= idx ; i++ )
    {
        cout << pr[i] << " ";
    }*/
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    build();
    return 0;
}