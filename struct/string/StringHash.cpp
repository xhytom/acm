using u64 = unsigned long long;
static const u64 P = (1ull << 61) - 1;
uniform_int_distribution<u64> dist(P / 2, P - 1);
const u64 base = dist(mrand);

struct HashValue {
    u64 x;
    HashValue(u64 _x) {
        x = _x;
    } 

    HashValue() {
        x = 0;
    }
    
    friend HashValue operator+(const HashValue &lhs, const HashValue &rhs) {
        u64 c = lhs.x + rhs.x;
        if (c >= P) {
            c -= P;
        }
        return {c};
    }

    friend HashValue operator-(const HashValue &lhs, const HashValue &rhs) {
        u64 c = lhs.x + P - rhs.x;
        if (c >= P) {
            c -= P;
        }
        return {c};
    }

    friend HashValue operator*(const HashValue &lhs, const HashValue &rhs) {
        __uint128_t c = __uint128_t(lhs.x) * rhs.x;
        return HashValue(c >> 61) + HashValue(c & P);
    }

    friend std::ostream &operator<<(std::ostream &os, HashValue value) {
        return os << value.x;
    }

    friend bool operator==(const HashValue &lhs, const HashValue &rhs) {
        return lhs.x == rhs.x;
    }
     friend bool operator!=(const HashValue &lhs, const HashValue &rhs) {
        return lhs.x != rhs.x;
    }
};



const struct HashBase {
    std::vector<HashValue> fact;
    HashBase(int n) : fact(n + 5){
        fact[0] = HashValue(1);
        for (int i = 1; i <= n + 3; i++) {
            fact[i] = fact[i - 1] * HashValue(base);
        }
    }
    HashBase(){}    
}B(N);

struct DynamicStringHash {
    std::vector<HashValue> pre, bac;
    DynamicStringHash() : pre(1){
        pre.reserve(N);
        pre[0] = {0};
    }
    int size() {
        return pre.size();
    }
    void push_back(char ch) {
        pre.push_back(pre.back() * HashValue(base) + HashValue(ch)) ;
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return (pre[r] - (pre[l - 1] * B.fact[r - l + 1]));
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[r - l + 1] + get(l, r));
    }
};

struct StringHash {
    std::vector<HashValue> pre, bac;
    StringHash(std::string s) : pre(s.size() + 1), bac(s.size() + 1){
        pre[0] = 0;
        for (int i = 1; i < s.size(); i++) {
            pre[i] = (pre[i - 1] * HashValue(base) + HashValue(s[i]));
        }
        bac[s.size()] = 0;
        for (int i = s.size() - 1; i >= 0; i--) {
            bac[i] = (bac[i + 1] * HashValue(base) + HashValue(s[i]));
        }

    }
    HashValue getrev(int l, int r) {
        if (l > r) return HashValue();
        return bac[l] - (bac[r + 1] * B.fact[r - l + 1]);
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return pre[r] - (pre[l - 1] * B.fact[r - l + 1]);
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[r - l + 1] + get(l, r));
    }
    HashValue concatrev(HashValue value, int l, int r) {
        return (value * B.fact[r - l + 1] + getrev(l, r));
    }
};