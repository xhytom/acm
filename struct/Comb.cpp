struct Comb {
    vector<ll> _fact, _invfact;
    Comb(int n) : _fact(n + 1), _invfact(n + 1) {
        _fact[0] = _invfact[0] = 1;
        for (int i = 1; i <= n; i++) {
            _fact[i] = (_fact[i - 1] * i) % MOD;
        }
        _invfact[n] = ksm(_fact[n], MOD - 2);
        for (int i = n - 1; i >= 1; i--) {
            _invfact[i] = _invfact[i + 1] * (i + 1) % MOD;
        }
    } 
    ll C(ll x, ll y) {
        return _fact[x] * _invfact[x - y] % MOD * _invfact[y] % MOD;
    }
    ll P(ll x, ll y) {
        return _fact[x] * _invfact[x - y] % MOD;
    }
    ll fact(ll x) {
        return _fact[x];
    }
    ll invfact(ll x) {
        return _invfact[x];
    }
    ll ksm(ll x, ll y) {
        ll ans = 1;
        while (y) {
            if (y & 1) ans = ans * x % MOD;
            y /= 2;
            x = x * x % MOD;
        }
        return ans;
    }
}Comb(N);