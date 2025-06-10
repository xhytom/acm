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
const ll N = 405;
//head
const ll MAXN = 1e6 + 5;
ll prime[MAXN], idxprime = 0;
bool isprime[MAXN];

void prime_build()
{
    for(int i = 2 ; i < MAXN ; i++ )
    {
        if(isprime[i] == 0)
        {
            prime[++idxprime] = i;
        }
        for(int j = 1 ; j <= idxprime && i * prime[j] < MAXN ; j++ )
        {
            isprime[i * prime[j]] = 1;
            if(i % prime[j] == 0) break;
        }
    }
    /*
    for(int i = 1 ; i <= idxprime ; i++ )
    {
        cout << prime[i] << " ";
    }
    */
}




int a[N];
int v[N], n1, n2;
int to[N], b[N];
int n;
vector<int> e[N];

bool find(int x)
{
    b[x] = true;
    for(auto y : e[x])
    {
        if(!v[y] || (!b[v[y]] && find(v[y])))
        {
            v[y] = x;
            return true;
        }
    }
    return false;
}

int match()
{
    int ans = 0;
    memset(v, 0 ,sizeof(v));
    for(int i = 1 ; i <= n1 ; i++ )
    {
        memset(b, 0, sizeof(b));
        if(find(i))
        {
            ++ans;
        }
    }
    return ans;
}

int main()
{  
    fastio
    //freopen("1.in","r",stdin);
    prime_build();
    cin >> n;
    for(int i = 1 ; i <= n ; i++ )
    {
        cin >> a[i];
        if(a[i] & 1)
        {
            to[i] = ++n1;
        }else{
            to[i] = ++n2;
        }
    }
    for(int i = 1 ; i <= n ; i++ )
    {
        for(int j = 1 ; j < i ; j++ )
        {
            if(isprime[a[i] + a[j]] == 0)
            {
                //cout << i <<":" << j << endl;
                if(a[i] & 1)
                {
                    e[to[i]].push_back(to[j]);
                }else{
                    e[to[j]].push_back(to[i]);
                }
            }
        }
    }
    cout << match() << "\n";
    return 0;
}