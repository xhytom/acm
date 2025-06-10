struct PAM {
    int sz, tot, last;
    vector<int> cnt, len, fail;
    vector<vector<int>> ch;
    vector<char> s;
    PAM(int n) : cnt(n + 5), ch(n + 5, vector<int>(30)), len(n + 5), fail(n + 5), s(n + 5) {
        clear();
    }
    int node(int l) {    // 建立一个新节点，长度为 l
        sz++;
        ch[sz].assign(30, 0);
        len[sz] = l;
        fail[sz] = cnt[sz] = 0;
        return sz;
    }
    void clear() {   // 初始化
        sz = -1;
        last = 0;
        s[tot = 0] = '$';
        node(0);
        node(-1);
        fail[0] = 1;
    }
    int getfail(int x) {   // 找后缀回文
        while (s[tot - len[x] - 1] != s[tot]) x = fail[x];
        return x;
    }
    void insert(char c)    // 建树
    {
        s[++tot] = c;
        int now = getfail(last);
        if (!ch[now][c - 'a'])
        {
            int x = node(len[now] + 2);
            fail[x] = ch[getfail(fail[now])][c - 'a'];
            ch[now][c - 'a'] = x;
        }
        last = ch[now][c - 'a'];
        cnt[last]++;
    }
};