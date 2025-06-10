std::vector<int> exkmp(string s) {
    std::vector<int> p(s.size());
    int n = s.size() - 1;
    int L = 1, R = 0;
    p[1] = 0;
    for(int i = 2 ; i <= n ; i++ ) {
        if(i > R) p[i] = 0;
        else p[i] = std::min(p[i - L + 1], R - i + 1);
        while (i + p[i] <= n && s[p[i] + 1] == s[i + p[i]]) ++p[i];
        if(i + p[i] - 1 > R) {
            L = i;
            R = i + p[i] - 1;
        }
    }
    return p;
}
