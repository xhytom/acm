using f32 = float;
using f64 = double;
using f128 = long double;
using a64 = double;
using a128 = long double;
using arc = double;

#define Vector Point
#define sp(x) cout << fixed << setprecision(x)

const f64 PI = acos(-1);
const f64 EPS = 1e-7;
const f64 INF = numeric_limits<f64>::max();

f64 fgcd(f64 a, f64 b) {
	return fabs(a) < EPS ? fabs(a) : fgcd(b, fmod(a, b));
}

template<class T, class S>
bool equal(T a, S b) {
	return -EPS < a - b && a - b < EPS;
}

template<class T>
int sign(T a) {
	if(-EPS < a && a < EPS) {
		return 0;
	}
	return a < 0 ? -1 : 1;
}

template<class T>
struct Point {
	T x, y;
	Point(T x_ = 0, T y_ = 0) : x(x_), y(y_) {}
	template<class U> operator Point<U>() {
		return Point<U> (U(x), U(y));
	}
	
	Point &operator += (Point p) & {return x += p.x, y += p.y, *this;}
	Point &operator += (T t) & {return x += t, y += t, *this;}
	Point &operator -= (Point p) & {return x -= p.x, y -= p.y, *this;}
	Point &operator -= (T t) & {return x -= t, y -= t, *this;}
	Point &operator *= (Point p) & {return x *= p.x, y *= p.y, *this;}
	Point &operator *= (T t) & {return x *= t, y *= t, *this;}
	Point &operator /= (T t) & {return x /= t, y /= t, *this;}
	Point operator - () const {return Point(-x, -y);}
	friend Point operator + (Point a, Point b) {return a += b;}
    friend Point operator + (Point a, T b) {return a += b;}
    friend Point operator - (Point a, Point b) {return a -= b;}
    friend Point operator - (Point a, T b) {return a -= b;}
    friend Point operator * (Point a, T b) {return a *= b;}
    friend Point operator * (T a, Point b) {return b *= a;}
    friend Point operator / (Point a, T b) {return a /= b;}
    friend T operator * (Point a, Point b) {return a.x * b.x + a.y * b.y;}
    friend T operator ^ (Point a, Point b) {return a.x * b.y - a.y * b.x;};
    
    friend bool operator < (Point a, Point b) {
    	return equal(a.x, b.x) ? a.y < b.y - EPS : a.x < b.x - EPS;
    }
    friend bool operator > (Point a, Point b) {return b < a;}
    friend bool operator == (Point a, Point b) {return !(a < b) && !(b < a);}
    friend bool operator != (Point a, Point b) {return a < b || b < a;}
    
    friend auto &operator>>(istream &is, Point &p) {
    	return is >> p.x >> p.y;
    }
    
    friend auto &operator<<(ostream &os, Point p) {
    	return os << "(" << p.x << ", " << p.y << ")";
    }
};


template<class T>
struct Line {
	Point<T> a, b;
	Line(Point<T> a_ = Point<T>(), Point<T> b_ = Point<T>()) : a(a_), b(b_) {}
	template<class U> operator Line<U>() {
		return Line<U>(Point<U>(a), Point<U>(b));
	}
	friend auto &operator << (ostream& os, Line l) {
		return os << "<" << l.a << ", " << l.b << ">";
	}
};

template<class T>
a128 atan(Point<T> p) { // 从 $x$ 负半轴逆时针排序 即 3 -> 4 -> 1 -> 2 左开右闭
	auto [x, y] = p;
	if(sign(x) < 0 && sign(y) == 0) {
		return 2 * PI;
	}
	if(sign(x) < 0 && sign(y) < 0) {
		return - atan2l(x, y) - PI / 2;
	}
	if(sign(x) == 0 && sign(y) < 0) {
		return PI / 2;
	}
	if(sign(x) > 0 && sign(y) < 0) {
		return PI * 3 / 2 - atan2l(x, y);
	}
	if(sign(x) >= 0 && sign(y) == 0) {
		return PI;
	}
	if(sign(x) > 0 && sign(y) > 0) {
		return PI * 3 / 2 - atan2l(x, y);
	}
	if(sign(x) == 0 && sign(y) > 0) {
		return 3 * PI / 2;
	}
	if(sign(x) < 0 && sign(y) > 0) {
		return 3 * PI / 2 - atan2l(x, y);
	}
	return 1e18;
}

template<class T>
a128 atanFromPosiX(Point<T> p) { // 从 $x$ 正半轴顺时针排序 即 4 -> 3 -> 2 -> 1 左开右闭
	auto [x, y] = p;
	if(sign(x) > 0 && sign(y) == 0) {
		return 0;
	}
	if(sign(x) > 0 && sign(y) < 0) {
		return atan2l(x, y) - PI / 2;
	}
	if(sign(x) == 0 && sign(y) < 0) {
		return PI / 2;
	}
	if(sign(x) < 0 && sign(y) < 0) {
		return PI * 3 / 2 - atan2l(x, y);
	}
	if(sign(x) < 0 && sign(y) == 0) {
		return PI;
	}
	if(sign(x) < 0 && sign(y) > 0) {
		return PI * 3 / 2 - atan2l(x, y);
	}
	if(sign(x) == 0 && sign(y) > 0) {
		return 3 * PI / 2;
	}
	if(sign(x) > 0 && sign(y) > 0) {
		return 3 * PI / 2 + atan2l(x, y);
	}
	return 1e18;
}

template<class T>
vector<Point<T>> sortByArgument(vector<Point<T>> vec) {
	sort(vec.begin() + 1, vec.end(), [&](Point<T> a, Point<T> b) {
		return atan(a) < atan(b);
	});
	return vec;
}


template<class T> 
T cross(Point<T> a, Point<T> b) {
    return a ^ b;
}

template<class T> 
T cross(Point<T> p1, Point<T> p2, Point<T> p0) { // p0 -> p1, p0 -> p2
    return (p1 - p0) ^ (p2 - p0);
}

template<class T> 
T dot(Point<T> a, Point<T> b) {
    return a * b;
}

template<class T> 
T dot(Point<T> p1, Point<T> p2, Point<T> p0) { // p0 -> p1, p0 -> p2
    return (p1 - p0) * (p2 - p0);
}

template <class T> 
T dis2(T x1, T y1, T x2, T y2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

template <class T> 
T dis2(Point<T> a, Point<T> b) {
    return dis2(a.x, a.y, b.x, b.y);
}

template <class T> 
f64 dis(T x1, T y1, T x2, T y2) {
    return sqrt(dis2(x1, y1, x2, y2));
}

template <class T> 
f64 dis(Point<T> a, Point<T> b) {
    return dis(a.x, a.y, b.x, b.y);
}

template<class T>
f64 length(Vector<T> v) {
	return sqrt(v.x * v.x + v.y * v.y);
}

Vector<f64> standardize(Vector<f64> v) {
	return v / length(v);
}

f64 toDeg(a64 x) { // 弧度转角度
    return x * 180 / PI;
}

a64 toArc(f64 x) { // 角度转弧度
    return PI / 180 * x;
}

a64 getArc(f64 a, f64 b, f64 c) {
	return acos((a * a + b * b - c * c) / (2.0 * a * b));	
}

f64 getDeg(f64 a, f64 b, f64 c) {
	return toDeg(getArc(a, b, c));	
}

template<class T>
a64 getArc(Point<T> a, Point<T> b) {
	return fabs(atan2(abs(a ^ b), a * b));
}

template<class T>
f64 getDeg(Point<T> a, Point<T> b) {
	return toDeg(getArc(a, b));
}

Point<f64> rotate(Point<f64> p, a64 rad) {
    return {p.x * cos(rad) - p.y * sin(rad), p.x * sin(rad) + p.y * cos(rad)};
}

template<class T> Point<T> rotate(Point<T> p, Point<T> base) { // p 绕 base 逆时针 90
    Vector<T> vec = p - base;
    return Point(-vec.y, vec.x);
}

Point<f64> rotate(Point<f64> p, Point<f64> base, a64 rad) {
    f64 x = (p.x - base.x) * cos(rad) + (p.y - base.y) * sin(rad) + base.x;
    f64 y = (base.x - p.x) * sin(rad) + (p.y - base.y) * cos(rad) + base.y;
    return {x, y};
}

template<class T> 
bool onLine(Point<T> a, Point<T> b, Point<T> c) {
    return sign(cross(b, a, c)) == 0;
}

template<class T> 
bool onLine(Point<T> p, Line<T> l) {
    return onLine(p, l.a, l.b);
}

template<class T> 
bool pointOnLineLeft(Point<T> p, Line<T> l) {
    return cross(l.b, p, l.a) > 0;
}

template<class T> 
bool pointOnLineSide(Point<T> p1, Point<T> p2, Line<T> vec) {
    T val = cross(p1, vec.a, vec.b) * cross(p2, vec.a, vec.b);
    return sign(val) == 1;
}

template<class T> 
bool pointNotOnLineSide(Point<T> p1, Point<T> p2, Line<T> vec) {
    T val = cross(p1, vec.a, vec.b) * cross(p2, vec.a, vec.b);
    return sign(val) == -1;
}

Point<f64> lineIntersection(Line<f64> l1, Line<f64> l2) {
	return l1.a + cross(l2.b, l1.a, l2.a) / cross(l2.b - l2.a, l1.a - l1.b) * (l1.b - l1.a);
}

template<class T> 
bool lineParallel(Line<T> p1, Line<T> p2) {
    return sign(cross(p1.a - p1.b, p2.a - p2.b)) == 0;
}
template<class T> 
bool lineVertical(Line<T> p1, Line<T> p2) {
    return sign(dot(p1.a - p1.b, p2.a - p2.b)) == 0;
}
template<class T> 
bool lineSame(Line<T> l1, Line<T> l2) {
    return lineParallel(Line{l1.a, l2.b}, {l1.b, l2.a}) &&
           lineParallel(Line{l1.a, l2.a}, {l1.b, l2.b}) && lineParallel(l1, l2);
}

f64 disToLine(Point<f64> p, Line<f64> l) {
    Point<f64> ans = lineIntersection({p, p + rotate(l.a, l.b)}, l);
    return dis(p, ans);
}

f64 dis2ToLine(Point<f64> p, Line<f64> l) {
    Point<f64> ans = lineIntersection({p, p + rotate(l.a, l.b)}, l);
    return dis2(p, ans);
}

template<class T>
Point<f64> nearestToLine(Point<T> p, Line<T> l) {
    Point<f64> ans = lineIntersection({p, p + rotate(l.a, l.b)}, l);
    return ans;
}

template<class T>
bool pointOnSegment(Point<T> p, Line<T> l) {
    return sign(cross(p, l.a, l.b)) == 0 && min(l.a.x, l.b.x) <= p.x && p.x <= max(l.a.x, l.b.x) &&
           min(l.a.y, l.b.y) <= p.y && p.y <= max(l.a.y, l.b.y);
}
template<class T> 
bool pointOnSegmentNonStrict(Point<T> p, Line<T> l) {
    return pointOnSegment(p, l) && min(l.a.x, l.b.x) < p.x && p.x < max(l.a.x, l.b.x) &&
           min(l.a.y, l.b.y) < p.y && p.y < max(l.a.y, l.b.y);
}

Point<f64> nearestToSegment(Point<f64> p, Line<f64> l) {
    if (sign(dot(p, l.b, l.a)) == -1) { // 特判到两端点的距离
        return l.a;
    } else if (sign(dot(p, l.a, l.b)) == -1) {
        return l.b;
    }
    return nearestToLine(p, l);
}

f64 disToSegment(Point<f64> p, Line<f64> l) {
    if (sign(dot(p, l.b, l.a)) == -1) { // 特判到两端点的距离
        return dis(p, l.a);
    } else if (sign(dot(p, l.a, l.b)) == -1) {
        return dis(p, l.b);
    }
    return disToLine(p, l);
}

Point<f64> project(Point<f64> p, Line<f64> l) { // 点在直线投影
    Vector<f64> v = l.b - l.a;
    f64 r = dot(v, p - l.a) / length(v);
    return l.a + v * r;
}

template<class T> 
Line<T> midSegment(Line<T> l) {
    Point<T> mid = (l.a + l.b) / 2;
    return {mid, mid + rotate(l.a, l.b)};
}

template<class T> 
tuple<int, Point<T>, Point<T>> segmentIntersection(Line<T> l1, Line<T> l2) {
    auto [s1, e1] = l1;
    auto [s2, e2] = l2;
    auto A = max(s1.x, e1.x), AA = min(s1.x, e1.x);
    auto B = max(s1.y, e1.y), BB = min(s1.y, e1.y);
    auto C = max(s2.x, e2.x), CC = min(s2.x, e2.x);
    auto D = max(s2.y, e2.y), DD = min(s2.y, e2.y);
    if (A < CC || C < AA || B < DD || D < BB) {
        return {0, {}, {}};
    }
    if (sign(cross(e1 - s1, e2 - s2)) == 0) { // parallel
        if (sign(cross(s2, e1, s1)) != 0) {
            return {0, {}, {}};
        }
        Point<T> p1(max(AA, CC), max(BB, DD));
        Point<T> p2(min(A, C), min(B, D));
        if (!pointOnSegment(p1, l1)) {
            swap(p1.y, p2.y);
        }
        if (p1 == p2) {
            return {3, p1, p2};
        } else {
            return {2, p1, p2};
        }
    } 
    auto cp1 = cross(s2 - s1, e2 - s1);
    auto cp2 = cross(s2 - e1, e2 - e1);
    auto cp3 = cross(s1 - s2, e1 - s2);
    auto cp4 = cross(s1 - e2, e1 - e2);
    if (sign(cp1 * cp2) == 1 || sign(cp3 * cp4) == 1) {
        return {0, {}, {}};
    }
    // 使用下方函数时请使用浮点数
    Point<f64> p = lineIntersection(l1, l2);
    if (sign(cp1) != 0 && sign(cp2) != 0 && sign(cp3) != 0 && sign(cp4) != 0) {
        return {1, p, p};
    } else {
        return {3, p, p};
    }
}

template<class T>
struct Circle {
	Point<T> o;
	T r;
	
	Circle(Point<T> o_ = Point<T>(), T r_ = 0) : o(o_), r(r_) {}
	template<class U> operator Circle<U>() {
		return Circle<U> (Point<U>(o), U(r));
	}	
};

pair<Point<f64>, f64> pointToCircle(Point<f64> p, Circle<f64> c) {
	Point<f64> U = c.o, V = c.o;
	f64 d = dis(p, c.o);
	if(sign(d) == 0) {
		return {c.o, 0.};
	}
	Vector<f64> v = standardize(p - c.o);
	U += v * c.r, V -= v * c.r;
	if(sign(dis(c.o, U) - dis(c.o, V)) == 1) {
		return {V, dis(c.o, V)};
	} else {
		return {U, dis(c.o, U)};
	}
}

Point<f64> radToPoint(Circle<f64> c, a64 rad) {
	Vector<f64> v = {c.r, 0};
	return c.o + rotate(v, rad);
}

tuple<int, Point<f64>, Point<f64>> circleIntersection(Circle<f64> c1, Circle<f64> c2) {
	f64 d = dis(c1.o, c2.o);
	if(sign(c1.r - c2.r) == 1) {
		swap(c1, c2);
	}
	if(sign(d - c1.r - c2.r) == 0) {
		return {1, c1.r + standardize(c2.o - c1.o) * c1.r, c1.r + standardize(c2.o - c1.o) * c1.r};
	} else if(sign(d - c1.r - c2.r) == 1) {
		return {0, {}, {}};
	} else if(sign(d + c1.r - c2.r) == -1) {
		return {0, {}, {}};
	} 
	Vector<f64> v = c2.o - c1.o;
	a64 init = atanFromPosiX(v);
	a64 arc = getArc(c1.r, d, c2.r);
	return {2, c1.o + rotate(Vector(c1.r, 0.), arc + init), c1.o + rotate(Vector(c1.r, 0.), init - arc)};
}