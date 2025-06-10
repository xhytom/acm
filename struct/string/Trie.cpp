struct Trie {
    std::vector<std::vector<int>> nxt;
    int idx = -1;

    Trie() {
        newnode();
        newnode();
    }

    int newnode() {
        nxt.push_back(std::vector<int>(26, 0));
        ++idx;
        return idx;
    }

    void insert(const string &s) {
        int now = 1;
        for (auto x : s) {
            if (nxt[now][x - 'a'] == 0) {
                nxt[now][x - 'a'] = newnode();
            }
            now = nxt[now][x - 'a'];
        }
    }

    std::string query(const string &s) {
        int now = 1;
        for (auto x : s) {
            if (nxt[now][x - 'a'] == 0) {
                return "";
            }
            now = nxt[now][x - 'a'];
        }
    }

};