void induced_sort(const vector<int> &vec, int val_range, vector<int> &SA, const vector<bool> &sl, const vector<int> &lms_idx) {
    vector<int> l(val_range, 0), r(val_range, 0);
    for (int c : vec) {
        if (c + 1 < val_range) ++l[c + 1];
        ++r[c];
    }
    partial_sum(l.begin(), l.end(), l.begin());
    partial_sum(r.begin(), r.end(), r.begin());
    fill(SA.begin(), SA.end(), -1);
    for (int i = lms_idx.size() - 1; i >= 0; --i)
        SA[--r[vec[lms_idx[i]]]] = lms_idx[i];
    for (int i : SA)
        if (i >= 1 && sl[i - 1]) {
            SA[l[vec[i - 1]]++] = i - 1;
        }
    fill(r.begin(), r.end(), 0);
    for (int c : vec)
        ++r[c];
    partial_sum(r.begin(), r.end(), r.begin());
    for (int k = SA.size() - 1, i = SA[k]; k >= 1; --k, i = SA[k])
        if (i >= 1 && !sl[i - 1]) {
            SA[--r[vec[i - 1]]] = i - 1;
        }
}
vector<int> SA_IS(const vector<int> &vec, int val_range) {
    const int n = vec.size();
    vector<int> SA(n), lms_idx;
    vector<bool> sl(n);
    sl[n - 1] = false;
    for (int i = n - 2; i >= 0; --i) {
        sl[i] = (vec[i] > vec[i + 1] || (vec[i] == vec[i + 1] && sl[i + 1]));
        if (sl[i] && !sl[i + 1]) lms_idx.push_back(i + 1);
    }
    reverse(lms_idx.begin(), lms_idx.end());
    induced_sort(vec, val_range, SA, sl, lms_idx);
    vector<int> new_lms_idx(lms_idx.size()), lms_vec(lms_idx.size());
    for (int i = 0, k = 0; i < n; ++i)
        if (!sl[SA[i]] && SA[i] >= 1 && sl[SA[i] - 1]) {
            new_lms_idx[k++] = SA[i];
        }
    int cur = 0;
    SA[n - 1] = cur;
    for (size_t k = 1; k < new_lms_idx.size(); ++k) {
        int i = new_lms_idx[k - 1], j = new_lms_idx[k];
        if (vec[i] != vec[j]) {
            SA[j] = ++cur;
            continue;
        }
        bool flag = false;
        for (int a = i + 1, b = j + 1;; ++a, ++b) {
            if (vec[a] != vec[b]) {
                flag = true;
                break;
            }
            if ((!sl[a] && sl[a - 1]) || (!sl[b] && sl[b - 1])) {
                flag = !((!sl[a] && sl[a - 1]) && (!sl[b] && sl[b - 1]));
                break;
            }
        }
        SA[j] = (flag ? ++cur : cur);
    }
    for (size_t i = 0; i < lms_idx.size(); ++i)
        lms_vec[i] = SA[lms_idx[i]];
    if (cur + 1 < (int)lms_idx.size()) {
        auto lms_SA = SA_IS(lms_vec, cur + 1);
        for (size_t i = 0; i < lms_idx.size(); ++i) {
            new_lms_idx[i] = lms_idx[lms_SA[i]];
        }
    }
    induced_sort(vec, val_range, SA, sl, new_lms_idx);
    return SA;
}
template <class T>
vector<int> suffix_array(const T &s, const int LIM = 128) {
    vector<int> vec(s.size() + 1);
    copy(begin(s), end(s), begin(vec));
    vec.back() = 0 ;
    // vec.back() = '$';
    auto ret = SA_IS(vec, LIM);
    ret.erase(ret.begin());
    return ret;
}
vector<int> getRank(const vector<int> &sa) {
    vector<int> rk(sa.size());
    for (int i = 0 ; i < sa.size(); i++) {
        rk[sa[i]] = i;
    }
    return rk;
}
template <class T>
vector<int> getHeight(const T &s, const vector<int> &sa) {
    int n = s.size(), k = 0;
    vector<int> ht(n), rank(n);
    for (int i = 0; i < n; i++) rank[sa[i]] = i;
    for (int i = 0; i < n; i++, k ? k-- : 0) {
        if (rank[i] == n - 1) {
            k = 0;
            continue;
        }
        int j = sa[rank[i] + 1];
        while (i + k < n && j + k < n && s[i + k] == s[j + k]) ++ k;
        ht[rank[i] + 1] = k;
    }
    ht[0] = 0;
    return ht;
}
template <class T>
vector<vector<int>> buildLCP(const T &s, const vector<int> ht) {
    vector<vector<int>> st;
    int n = s.size() - 1;
    int LOG = __lg(n) + 1;
    st.resize(LOG);
    st[0].resize(n + 1); 
    for(int i = 1 ; i <= n ; i++ )
    {
        st[0][i] = ht[i];
    }
    for(int j = 1 ; j <= LOG ; j++ )
    {
        st[j].resize(n + 1);
        for(int i = 1 ; i + (1 << j) - 1 <= n ; i++ )
        {
            st[j][i] = min(st[j - 1][i], st[j - 1][i + (1ll << j - 1)]);
        }
    }
    return st;
}
void use() {
    vector<vector<int>> st;
    vector<int> rk;
    int n;
    int u, v;
    function<int(int, int)> lcp = [&](int u, int v)
    {
        if(u == v) return n - u + 1;
        if(rk[u] > rk[v]) swap(u, v);
        int l = rk[u] + 1, r = rk[v];
        int len = __lg(r - l + 1);
        return min(st[len][l], st[len][r - (1 << len) + 1]);
    };
}