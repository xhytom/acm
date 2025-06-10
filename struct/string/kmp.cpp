std::vector<int> kmp(const string &s) {
    std::vector<int> nxt(s.size());
    int n = s.size() - 1;
    nxt[1] = 0;
    int j = 0;
    for (int i = 2; i <= n; i++) {
        while (s[j + 1] != s[i] && j != 0) j = nxt[j];
        if (s[j + 1] == s[i]) j++;
        nxt[i] = j;
    }
    return nxt;
}