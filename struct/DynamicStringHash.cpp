struct HashValue {
    int x = 0, y = 0;
    HashValue(int _x) {
        x = y = _x;
    } 
    HashValue(int _x, int _y) {
        x = _x;
        y = _y;
    }
    HashValue(){}
    
    friend HashValue operator+(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x + rhs.x, lhs.y + rhs.y};
    }

    friend HashValue operator%(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x % rhs.x, lhs.y % rhs.y};
    }

    friend HashValue operator*(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x * rhs.x, lhs.y * rhs.y};
    }

    friend HashValue operator-(const HashValue &lhs, const HashValue &rhs) {
        return {lhs.x - rhs.x, lhs.y - rhs.y};
    }

    friend std::ostream &operator<<(std::ostream &os, HashValue value) {
        return os << value.x << " " << value.y;
    }

    friend bool operator==(const HashValue &lhs, const HashValue &rhs) {
        return lhs.x == rhs.x && lhs.y == rhs.y;
    }
};



const struct HashBase {
    std::vector<HashValue> fact;
    HashValue P, base;
    int nextPrime(int n) {
        int flag = 1;
        while (flag) {
            n++;
            flag = 0;
            for (int i = 2; i * i <= n; i++) {
                if (n % i == 0) {
                    flag = 1;
                    break;
                }
            }
        }
        return n;
    }
    HashBase(int n) : fact(n + 5){
        const int P1 = nextPrime(rnd(900000000, 1000000000));
        const int P2 = nextPrime(rnd(900000000, 1000000000));
        P = HashValue(P1, P2);
        const int base1 = nextPrime(rnd(200, 10000));
        const int base2 = nextPrime(rnd(200, 10000));
        base = HashValue(base1, base2);
        fact[0] = HashValue(1);
        for (int i = 1; i <= n + 3; i++) {
            fact[i] = fact[i - 1] * base % P;
        }
    }
    HashBase(){}    
}B(N);

struct DynamicStringHash {
    std::vector<HashValue> pre, bac;
    DynamicStringHash() : pre(1){
        pre.reserve(N);
        pre[0] = {0, 0};
    }
    int size() {
        return pre.size();
    }
    void push_back(char ch) {
        pre.push_back(pre.back() * B.base % B.P + HashValue(ch) % B.P) ;
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return ((pre[r] - (pre[l - 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + get(l, r)) % B.P;
    }
};

struct StringHash {
    std::vector<HashValue> pre, bac;
    StringHash(std::string s) : pre(s.size() + 1), bac(s.size() + 1){
        pre[0] = {0, 0};
        for (int i = 1; i < s.size(); i++) {
            pre[i] = (pre[i - 1] * B.base % B.P + HashValue(s[i])) % B.P;
        }
        bac[s.size()] = {0, 0};
        for (int i = s.size() - 1; i >= 0; i--) {
            bac[i] = (bac[i + 1] * B.base % B.P + HashValue(s[i])) % B.P;
        }

    }
    HashValue getrev(int l, int r) {
        if (l > r) return HashValue();
        return ((bac[l] - (bac[r + 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue get(int l, int r) {
        if (l > r) return HashValue();
        return ((pre[r] - (pre[l - 1] * B.fact[r - l + 1]) % B.P) + B.P) % B.P;
    }
    HashValue concat(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + get(l, r)) % B.P;
    }
    HashValue concatrev(HashValue value, int l, int r) {
        return (value * B.fact[std::abs(r - l) + 1] % B.P + getrev(l, r)) % B.P;
    }
};