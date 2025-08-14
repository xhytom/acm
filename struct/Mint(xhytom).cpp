using u64 = unsigned long long;
using u32 = unsigned;
using i64 = long long;
using i32 = signed;

template <i32 MOD>
u32 down(u32 x) { return x >= MOD ? x - MOD : x; }

template <i32 MOD>
struct MInt {
    u32 x;
    MInt () : x(0) {}
    MInt (u32 x) : x(x) {}
    friend MInt operator+(MInt a, MInt b) { return down<MOD>(a.x + b.x); }
    friend MInt operator-(MInt a, MInt b) { return down<MOD>(a.x - b.x + MOD); }
    friend MInt operator*(MInt a, MInt b) { return 1ULL * a.x * b.x % MOD; }
    friend MInt operator/(MInt a, MInt b) { return a * ~b; }
    friend MInt operator^(MInt a, i64 b) {
        MInt ans = 1;
        while (b) {
            if (b & 1) ans = ans * a;
            a = a * a;
            b /= 2;
        }
        return ans;
    }
    friend MInt operator~(MInt a) { return a ^ (MOD - 2); }
    // friend MInt operator~(MInt a) { return MInt(INV::inv(a.x)); }
    friend std::istream &operator>>(std::istream &in, MInt &a) { return in >> a.x; }
    friend std::ostream &operator<<(std::ostream &out, MInt a) { return out << a.x; }
    friend MInt operator-(MInt a) { return down<MOD>(MOD - a.x); }
    friend MInt &operator+=(MInt &a, MInt b) { return a = a + b; }
    friend MInt &operator-=(MInt &a, MInt b) { return a = a - b; }
    friend MInt &operator*=(MInt &a, MInt b) { return a = a * b; }
    friend MInt &operator/=(MInt &a, MInt b) { return a = a / b; }
    friend MInt &operator^=(MInt &a, long long b) { return a = a ^ b; }
    friend bool operator==(MInt a, MInt b) { return a.x == b.x; }
    friend bool operator!=(MInt a, MInt b) { return !(a == b); }
    friend bool operator<(MInt a, MInt b) { return a.x < b.x; }
};

const i32 MOD = 998244353;
using Z = MInt<MOD>;