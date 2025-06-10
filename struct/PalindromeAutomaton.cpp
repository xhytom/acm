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