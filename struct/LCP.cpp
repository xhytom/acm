using i64 = int64_t;
using u64 = uint64_t;
struct LongestCommonPrefix {
    int n;
    vector<int> p, rank;
    vector<vector<int>> st;
    LongestCommonPrefix(const string &s) : n(s.size()), p(n), rank(n) {
        int k = 0;
        vector<int> q, count;

        for (int i = 0; i < n; i += 1)
            p[i] = i;

        sort(p.begin(), p.end(), [&](int i, int j) {
            return s[i] < s[j];
        });

        for (int i = 0; i < n; i += 1)
            rank[p[i]] = i and s[p[i]] == s[p[i - 1]] ? rank[p[i - 1]] : k++;

        for (int m = 1; m < n; m *= 2) {
            q.resize(m);

            for (int i = 0; i < m; i += 1)
                q[i] = n - m + i;

            for (int i : p)
                if (i >= m)
                    q.push_back(i - m);

            count.assign(k, 0);

            for (int i : rank)
                count[i] += 1;

            for (int i = 1; i < k; i += 1)
                count[i] += count[i - 1];

            for (int i = n - 1; i >= 0; i -= 1)
                p[count[rank[q[i]]] -= 1] = q[i];

            auto cur = rank;
            cur.resize(2 * n, -1);
            k = 0;

            for (int i = 0; i < n; i += 1)
                rank[p[i]] = i and cur[p[i]] == cur[p[i - 1]] and
                             cur[p[i] + m] == cur[p[i - 1] + m]
                             ? rank[p[i - 1]]
                             : k++;
        }

        st.emplace_back(n);

        for (int i = 0, k = 0; i < n; i += 1) {
            if (not rank[i])
                continue;

            k = max(k - 1, 0LL);
            int j = p[rank[i] - 1];

            while (i + k < n and j + k < n and s[i + k] == s[j + k])
                k += 1;

            st[0][rank[i]] = k;
        }

        for (int i = 1; (1 << i) < n; i += 1) {
            st.emplace_back(n - (1 << i) + 1);

            for (int j = 0; j <= n - (1 << i); j += 1)
                st[i][j] = min(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
        }
    }
    int get(int i, int j) {
        if (i == j)
            return n - i;

        if (i == n or j == n)
            return 0;

        i = rank[i];
        j = rank[j];

        if (i > j)
            swap(i, j);

        int k = 64 - __builtin_clzll(u64(j - i)) - 1;
        return min(st[k][i + 1], st[k][j - (1 << k) + 1]);
    }
};