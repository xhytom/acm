struct ACautomaton {
    vector<vector<int>> nxt;
    vector<int> fail;
    int idx = 0;
    ACautomaton() : nxt(1, vector<int>(26, 0)), fail(1){}

    int newnode() {
        int cur = ++idx;
        nxt.push_back(vector<int>(26, 0));
        fail.emplace_back(0);
        return cur;
    }
    void insert(const string &s) {
        int now = 0;
        for (auto x : s) {
            if (!nxt[now][x - 'a']) {
                nxt[now][x - 'a'] = newnode();
            }
            now = nxt[now][x - 'a'];
        }
    }
    void buildfail() {
        queue<int> q;
        for (int i = 0; i <= 25; i++) {
            if (nxt[0][i]) {
                fail[nxt[0][i]] = 0;
                q.push(nxt[0][i]);
            }
        }
        while (!q.empty()) {
            int now = q.front();
            q.pop();
            for (int i = 0; i <= 25; i++) {
                if (nxt[now][i]) {
                    fail[nxt[now][i]] = nxt[fail[now]][i];
                    q.push(nxt[now][i]);
                } else {
                    nxt[now][i] = nxt[fail[now]][i];
                }
            }
        }
    }
};// 推荐root = 0，root = 1会造成很多麻烦