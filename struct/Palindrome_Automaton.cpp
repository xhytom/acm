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
// #define int long long
#define pb push_back
#define eb emplace_back

const ll P1 = 999971, base1 = 101;
const ll P2 = 999973, base2 = 103;
// const ll N = 200005;
//head

const int N = 5e5 + 10, Sigma = 26;
char s[N];
int lastans, n;
struct Palindrome_Automaton {
    int ch[N][Sigma], fail[N], len[N], sum[N], cnt, last;
    Palindrome_Automaton() {
        cnt = 1;
        fail[0] = 1, fail[1] = 1, len[1] = -1;
    }
    int getfail(int x, int i) {
        while(i - len[x] - 1 < 0 || s[i - len[x] - 1] != s[i]) x = fail[x];
        return x;
    }
    void insert(char c, int i) {
        int x = getfail(last, i), w = c - 'a';
        if(!ch[x][w]) {
            len[++cnt] = len[x] + 2;
            int tmp = getfail(fail[x], i);
            fail[cnt] = ch[tmp][w];
            sum[cnt] = sum[fail[cnt]] + 1;
            ch[x][w] = cnt;
        }
        last = ch[x][w];
    }
    
} PAM;


signed main()
{
#ifdef localfreopen
    // freopen("1.in","r",stdin);
#endif
    fastio
    scanf("%s", s + 1);
    int last = 0;
    int n = strlen(s + 1);
    for(int i = 1 ; i <= n ; i++ )
    {

        s[i] = (s[i] - 97 + last) % 26 + 97;
        PAM.insert(s[i], i);
        int ans = PAM.sum[PAM.last];
        cout << ans << " ";
        last = ans;
    }
    return 0;
}
