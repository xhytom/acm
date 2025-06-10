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
#define rep(i,a,n) for(int i=a;i<n;i++)
#define per(i,a,n) for(int i=n-1;i>=a;i--)
#define fastio ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define multi int _;cin>>_;while(_--)
#define debug(x) cerr << #x << " = " << (x) << endl;
#define int long long
#define pb push_back
#define eb emplace_back
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
void test() {cerr << "\n";}
template<typename T, typename... Args> 
void test(T x, Args... args) {cerr << x << " ";test(args...);}
// const ll MOD = 998244353;
const ll MOD = 1e9+7;
int ksm(int x,int y){int ans=1;x%=MOD;while(y){if(y&1)ans=ans*x%MOD;x=x*x%MOD,y/=2;}return ans;}

const ll P1 = 999971, base1 = 101;
const ll P2 = 999973, base2 = 103;
const ll N = 200005;
//head

vector<int> kmp(string s)
{//string的形式为'#' + t1 + '#' + s
	int n = s.size() - 1;
	vector<int> nxt(s.size());
	int j = 0;
	for(int i = 2 ; i <= n ; i++ ){
		while(j && s[j + 1] != s[i]) j = nxt[j];
		if(s[j + 1] == s[i]) j++;
		nxt[i] = j;
	}
	return nxt;
}

signed main()
{  
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    string ss, t;
    cin >> ss >> t;
    string s = '#' + t + '#' + ss;
    auto nxt = kmp(s);
    for(int i = 1 ; i < s.size() ; i++ )
    {
    	if(nxt[i] == t.size())
    	{
    		cout << i - t.size() * 2 << "\n";
    	}
    }
    for(int i = 1 ; i <= t.size() ; i++ )
    {
    	cout << nxt[i] << " \n"[i == t.size()];
    }
    return 0;
}