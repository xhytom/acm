#include <bits/stdc++.h>
using namespace std;

#define int long long

struct P {
    int x, y;
    P(){}
    P(int x_, int y_) : x(x_), y(y_){}
    void read() {
        cin >> x >> y;
    }
};

void solve() {
    int a, b;
    cin >> a >> b;

    vector<P> ps(10);
    for(int i = 1; i <= 3; i++) {
        ps[i].read();
    }
    int dis = (ps[3].x + ps[3].y - ps[1].x - ps[1].y);
    double t = 1. * dis / a;

    if(a >= b) {
        cout << t << '\n';
        return;
    }

    if(ps[1].x > ps[3].x) {
        ps[1].x = 2 * ps[3].x - ps[1].x;
        ps[2].x = 2 * ps[3].x - ps[2].x;
    }

    if(ps[1].y > ps[3].y) {
        ps[1].y = 2 * ps[3].y - ps[1].y;
        ps[2].y = 2 * ps[3].y - ps[2].y;
    }

    if(ps[2].x < ps[1].x) {
        if(ps[2].y < ps[1].y) {
            double zj = 1. * (ps[2].x - ps[1].x + ps[2].y - ps[1].y) / (b - a);
            t = min(t, (dis - zj * a) / b + zj);
        } else if(ps[2].y < ps[3].y) {
            double tb = 1. * abs(ps[1].x - ps[2].x) / b;
            double ta = 1. * abs(ps[1].y - ps[2].y) / a;
            if(ta < tb) {
                double dt = 1. * (abs(ps[1].x - ps[2].x) - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[1].y - ps[2].y) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } else {
            int d = abs(ps[2].x - ps[1].x) + abs(ps[2].y - ps[3].y);
            double tb = 1. * d / b;
            double ta = 1. * abs(ps[3].y - ps[1].y) / a;
            if(ta < tb) {
                double dt = 1. * (d - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[3].y - ps[1].y) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } 
    } else if(ps[2].x < ps[3].x) {
        if(ps[2].y < ps[1].y) {
            double tb = 1. * abs(ps[1].y - ps[2].y) / b;
            double ta = 1. * abs(ps[1].x - ps[2].x) / a;
            if(ta < tb) {
                double dt = 1. * (abs(ps[1].y - ps[2].y) - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[1].x - ps[2].x) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } else if(ps[2].y < ps[3].y) {
            double xy = 1. * (ps[2].x - ps[1].x + ps[2].y - ps[1].y) / (b + a);
            t = min(t, (dis - xy * a) / b + xy);
        } else {
            double tb = 1. * abs(ps[2].y - ps[3].y) / b;
            double ta = 1. * (abs(ps[1].y - ps[3].y) + abs(ps[1].x - ps[2].x)) / a;
            if(ta < tb) {
                double dt = 1. * (abs(ps[2].y - ps[3].y) - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[1].y - ps[3].y) + abs(ps[1].x - ps[2].x) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } 
    } else {
        if(ps[2].y < ps[1].y) {
            int d = abs(ps[2].x - ps[3].x) + abs(ps[2].y - ps[1].y);
            double tb = 1. * d / b;
            double ta = 1. * abs(ps[1].x - ps[3].x) / a;
            if(ta < tb) {
                double dt = 1. * (d - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[3].x - ps[1].x) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } else if(ps[2].y < ps[3].y) {
            double tb = 1. * abs(ps[2].x - ps[3].x) / b;
            double ta = 1. * (abs(ps[1].y - ps[2].y) + abs(ps[1].x - ps[3].x)) / a;
            if(ta < tb) {
                double dt = 1. * (abs(ps[2].x - ps[3].x) - ta * b) / (b - a);
                t = min(t, (dis - (ta + dt) * a) / b + ta + dt);
            } else {
                double dt = 1. * (abs(ps[1].y - ps[2].y) + abs(ps[1].x - ps[3].x) - tb * a) / (b + a);
                t = min(t, (dis - (tb + dt) * a) / b + tb + dt);
            }
        } else {
            int d = abs(ps[2].x - ps[1].x) + abs(ps[2].y - ps[1].y);
            double xy = 1. * d / (a + b);
            t = min(t, (dis - xy * a) / b + xy);
        } 
    }
    cout << fixed << setprecision(10);
    cout << t << '\n';

}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    int t = 1;
    cin >> t;
    while(t--) {
        solve();
    }
}
