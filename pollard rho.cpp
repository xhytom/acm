using i64 = long long;
using i128 = __int128;
i64 power(i64 a, i64 b, i64 m) {
    i64 res = 1;
    for (; b; b >>= 1, a = i128(a) * a % m) {
        if (b & 1) {
            res = i128(res) * a % m;
        }
    }
    return res;
}
 
bool isprime(i64 p) {
    if (p < 2) {
        return 0;
    }
    i64 d = p - 1, r = 0;
    while (!(d & 1)) {
        r++;
        d >>= 1;
    }
    int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    for (auto a : prime) {
        if (p == a) {
            return true;
        }
        i64 x = power(a, d, p);
        if (x == 1 || x == p - 1) {
            continue;
        }
        for (int i = 0; i < r - 1; i++) {
            x = i128(x) * x % p;
            if (x == p - 1) {
                break;
            }
        }
        if (x != p - 1) {
            return false;
        }
    }
    return true;
}
 
mt19937 rng((unsigned int) chrono::steady_clock::now().time_since_epoch().count());
 
i64 pollard_rho(i64 x) {
    i64 s = 0, t = 0;
    i64 c = i64(rng()) % (x - 1) + 1;
    i64 val = 1;
    for (int goal = 1; ; goal <<= 1, s = t, val = 1) {
        for (int step = 1; step <= goal; step++) {
            t = (i128(t) * t + c) % x;
            val = i128(val) * abs(t - s) % x;
            if (step % 127 == 0) {
                i64 g = gcd(val, x);
                if (g > 1) {
                    return g;
                }
            }
        }
        i64 g = gcd(val, x);
        if (g > 1) {
            return g;
        }
    }
}

unordered_map<i64, int> getprimes(i64 x) {
    unordered_map<i64, int> p;
    function<void(i64)> get = [&](i64 x) {
        if (x < 2) {
            return;
        }
        if (isprime(x)) {
            p[x]++;
            return;
        }
        i64 mx = pollard_rho(x);
        get(x / mx);
        get(mx);
    };
    get(x);
    return p;
}
