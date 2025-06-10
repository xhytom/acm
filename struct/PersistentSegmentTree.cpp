struct Node {
    int v = 0;
    Node *ls = nullptr;
    Node *rs = nullptr;
};
  
void add(Node *pre, Node *&cur, int L, int R, int x) {
    cur = new Node();
    if (pre != nullptr) {
        cur->v = pre->v;
        cur->ls = pre->ls;
        cur->rs = pre->rs;
    }
    cur->v++;
    // test(L, R);
    if (L == R) {
        return;
    }
      
    int mid = (L + R) / 2;
    if (x <= mid) {
        add(pre == nullptr ? nullptr : pre->ls, cur->ls, L, mid, x);
    } else {
        add(pre == nullptr ? nullptr : pre->rs, cur->rs, mid + 1, R, x);
    }
}
  
  
int query(Node *pre, Node *bac, int L, int R, int p) {
    if (L == R) {
        return L;
    }
      
    int sum = ((bac == nullptr || bac->ls == nullptr) ? 0 : bac->ls->v)
            - ((pre == nullptr || pre->ls == nullptr) ? 0 : pre->ls->v);
              
    int mid = (L + R) / 2;
    if (sum < p) {
        return query((pre == nullptr ? nullptr : pre->rs), (bac == nullptr ? nullptr : bac->rs), mid + 1, R, p - sum);
    } else {
        return query((pre == nullptr ? nullptr : pre->ls), (bac == nullptr ? nullptr : bac->ls), L, mid, p);
    }
      
}
  
int querySize(Node *pre, Node *bac, int L, int R, int l, int r) {
    if (R < l || r < L) return 0;
    if (l <= L && R <= r) {
        return (bac == nullptr ? 0 : bac->v) - (pre == nullptr ? 0 : pre->v);
    }
      
    int mid = (L + R) / 2;
    int ans = 0;
    ans += querySize((pre == nullptr ? nullptr : pre->ls), (bac == nullptr ? nullptr : bac->ls), L, mid, l, r);
    ans += querySize((pre == nullptr ? nullptr : pre->rs), (bac == nullptr ? nullptr : bac->rs), mid + 1, R, l, r);
    return ans;
}