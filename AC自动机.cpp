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
ll gcd(ll a,ll b){ return b?gcd(b,a%b):a;}
mt19937 mrand(random_device{}());
int rnd(int x){ return mrand() % x; }
void test() {cerr << "\n";}
template<typename T, typename... Args> 
void test(T x, Args... args) {cerr << x << " ";test(args...);}
const ll MOD = 998244353;
//const ll MOD = 1e9+7;
const ll P1 = 9999971, base1 = 101;
const ll P2 = 9999973, base2 = 103;
const ll N = 200005;
//head

int n, idx;
struct Node{
    int fail, nxt[26], end;
}trie[150000];


string ss[155];
int cnt[155];

void add_string(string s, int num)
{
    int p = 0;
    for(int i = 0 ; i < s.size() ; i++ )
    {
        int x = s[i] - 'a';
        if(!trie[p].nxt[x])
        {
            trie[p].nxt[x] = ++idx;
        }
        p = trie[p].nxt[x];
    }
    trie[p].end = num;
}

void get_fail()
{
    queue<int> q;
    rep(i, 0, 26)
    {
        if(trie[0].nxt[i])
        {
            trie[trie[0].nxt[i]].fail = 0;
            q.push(trie[0].nxt[i]);
        }
    }
    while(!q.empty())
    {
        int x = q.front();
        q.pop();
        rep(i, 0, 26)
        {
            if(trie[x].nxt[i])
            {
                trie[trie[x].nxt[i]].fail = trie[trie[x].fail].nxt[i];
                q.push(trie[x].nxt[i]);
            }else{
                trie[x].nxt[i] = trie[trie[x].fail].nxt[i];
            }
        }
    }
}
void query_string(string s)
{
    int p = 0;
    for(int i = 0 ; i < s.size() ; i++ )
    {
        int x = s[i] - 'a';
        if(trie[p].nxt[x])
        {
            p = trie[p].nxt[x];
        }else{
            p = trie[trie[p].fail].nxt[x];
        }
        for(int i = p ; i ; i = trie[i].fail)
        {
            cnt[trie[i].end]++;
        }
        //cout << p << " \n"[i == s.size() - 1];
    }   
}

signed main()
{  
    fastio
    //freopen("1.in","r",stdin);
    string s;
    while(cin >> n)
    {
        if(n == 0) break;
        idx = 0;
        memset(trie, 0, sizeof(trie));
        memset(cnt, 0, sizeof(cnt));
        rep(i, 1, n + 1)
        {
            cin >> ss[i];
            add_string(ss[i], i);
        }
        get_fail();
        cin >> s;
        query_string(s);
        ll ans = *max_element(cnt + 1, cnt + n + 1); 
        cout << ans << endl;
        /*rep(i, 1, n + 1)
        {
            cout << cnt[i] << " \n"[i == n];
        }*/
        rep(i, 1, n + 1)
        {
            if(cnt[i] == ans) cout << ss[i] << endl;
        }
        /*for(int i = 1 ; i <= idx ; i++ )
        {
            cout << i << "->" << trie[i].fail << "\n"; 
        }*/
    }
    
    return 0;
}