//
// Created by Armeria on 2020/3/14.
// UNFINISHED
//

#include <cstdio>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 100 + 5;
const double eps = 1e-12;

#define X_AXIS 0
#define Y_AXIS 1

int fsign(const double x) { return fabs(x) < eps ? 0 : (x < 0 ? -1 : 1); }

bool equals(const double a, const double b) { return fsign(a - b) == 0; }

struct Point {
    double x, y;

    Point() : x(0), y(0) {}

    Point(double x, double y) : x(x), y(y) {}

    bool operator <(const Point &a) const {
        return equals(x, a.x) ? y < a.y : x < a.x;
    }

    bool operator ==(const Point &a) {
        return equals(x, a.x) && equals(y, a.y);
    }

    bool operator !=(const Point &a) {
        return !equals(x, a.x) || !equals(y, a.y);
    }

} points[N];

typedef Point Vector;

Vector operator -(const Vector v) {
    return Vector(-v.x, -v.y);
}

Vector operator +(const Vector a, const Vector b) {
    return Vector(a.x + b.x, a.y + b.y);
}

Vector operator -(Vector a, Vector b) {
    return Vector(a.x - b.x, a.y - b.y);
}

Vector operator *(const Vector a, const double p) {
    return Vector(a.x * p, a.y * p);
}

Vector operator /(const Vector a, const double p) {
    return Vector(a.x / p, a.y / p);
}

double dot(Vector a, Vector b) { return a.x * b.x + a.y * b.y; }

double cross(Vector a, Vector b) { return a.x * b.y - a.y * b.x; }

double length(Vector a) { return sqrt(dot(a, a)); }

Vector rotate(Vector a, double rad) {
    return Vector(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad));
}

Vector normal(Vector a) {
    int l = length(a);
    return Vector(-a.y / l, a.x / l);
}

Point getLineIntersection(Point p, Vector v, Point q, Vector u) {
    Vector w = p - q;
    double t = cross(u, w) / cross(v, u);
    return p + v * t;
}

Point getLineIntersectionPts(Point a, Point b, Point c, Point d) {
    return getLineIntersection(a, b - a, c, d - c);
}

double getParameterT(Point p, Point a, Vector u) {
    return (p.x - a.x) / u.x;
}

double distance2Line(Point p, Point a, Point b) {
    Vector v1 = b - a, v2 = p - a;
    return fabs(cross(v1, v2) / length(v1));
}

double angle(Vector a, Vector b) {
    return acos(dot(a, b) / length(a) / length(b));
}

double getPerpendicularBisectorIntersection(Point a, Point b, int ParallelTo, double p) {
    Vector n0 = normal(a - b);
    if (ParallelTo == X_AXIS) {
        return getLineIntersection(Point(0, p), Vector(1, 0), (a + b) / 2.0, n0).x;
    }
    else {  // ParallelTo == Y_AXIS
        return getLineIntersection(Point(p, 0), Vector(0, 1), (a + b) / 2.0, n0).x;
    }
}

void partition(pair<double, pair<int, int>> pairs[], double xMid, int &i, int &j, int n) {
    i = 0;
    j = n - 1;
    while (i <= j) {
        while (pairs[i].first < xMid) i++;
        while (pairs[j].first > xMid) j--;
        if (i < j) {
            swap(pairs[i], pairs[j]);
            i++;
            j--;
        }
    }
}

double getKthSplit(pair<double, pair<int, int>> pairs[], int n, int k) {
    static const int c = 5;

    if (n < c) {
        sort(pairs, pairs + n);
        return pairs[k].first;
    }

    int s = 0;
    for (int i = 0; i < n - c + 1; i += c) {
        sort(pairs + i, pairs + i + c);
        swap(pairs[s++], pairs[i + (c >> 1)]);
    }

    double xMid = getKthSplit(pairs, s, ((s - 1) >> 1));
    pair<double, pair<int, int>> cpy[N];
    int m = 0, midNum = 0;

    int i, j;
    partition(pairs, xMid, i, j, n);
    for (int p = 0; p < n; p++) {
        if (pairs[p].first != xMid)
            cpy[m++] = pairs[p];
        else
            midNum++;
    }
    partition(cpy, xMid, i, j, m);


    int leftNum = i, rightNum = m - i;
    if (leftNum > k) return getKthSplit(cpy, i, k);
    else if (leftNum + midNum > k) return xMid;
    else return getKthSplit(cpy + leftNum, rightNum, k - midNum - leftNum);
}

void chooseNearerPoint(const pair<Point, Point> pair, const Point ref, Point filter[], int &s) {
    Point chosen;
    if (length(pair.first - ref) > length(pair.second - ref))
        chosen = pair.first;
    else
        chosen = pair.second;
    if (s == 0 || (s > 0 && chosen != filter[0])) {
        filter[s++] = chosen;
    }

}

pair<double, double> constraintYCoordianteCenter(Point points[], int n, double y) {
    if (n == 1) return make_pair(points[0].x, abs(points[0].y - y));
    else if (n == 2) {
        if (equals(points[0].x, points[1].x)) {
            return make_pair(points[0].x, max(fabs(points[0].y - y), fabs(points[1].y - y)));
        }
        double x = getPerpendicularBisectorIntersection(points[0], points[1], X_AXIS, y);
        Point c = Point(x, y);
        return make_pair(x, max(length(points[0] - c), length(points[1] - c)));
    }

    int pairCount = 0;
    static pair<double, pair<int, int>> pairs[N];
    for (int i = 0; i < n; i += 2) {
        int second = i + 1 == n ? 0 : i + 1;
        pair<Point, Point> cur = make_pair(points[i], points[second]);
        pairs[pairCount++] =
                make_pair(getPerpendicularBisectorIntersection(cur.first, cur.second, X_AXIS, y),
                          make_pair(i, second));
    }

    double xMid;
    if (pairCount & 1) {
        xMid = getKthSplit(pairs, pairCount, pairCount >> 1);
    }
    else {
        double mid0 = getKthSplit(pairs, pairCount, (pairCount - 1) >> 1);
        double mid1 = getKthSplit(pairs, pairCount, pairCount >> 1);
        xMid = (mid0 + mid1) / 2;
    }
    double farDistance = 0;
    Point pointMid = Point(xMid, y), farPoint;
    for (int i = 0; i < n; i++) {
        double len = length(points[i] - pointMid);
        if (farDistance < len) {
            farDistance = len;
            farPoint = points[i];
        }
    }

    static Point filter[N];
    int s = 0;
    if (farPoint.x > xMid) {
        int i;
        bool isDiscarded[N];
        for (i = 0; i < pairCount && pairs[i].first < xMid; i++) {
            pair<Point, Point> pair = make_pair(points[pairs[i].second.first], points[pairs[i].second.second]);
            if (length(pair.first - pointMid) < length(pair.second - pointMid)) {
                isDiscarded[pairs[i].second.first] = true;
            }
            else {
                isDiscarded[pairs[i].second.second] = true;
            }
        }
        for (int i = 0; i < n; i++) if (!isDiscarded[i]) filter[s++] = points[i];
    }
    else {
        int i;
        bool isDiscarded[N];
        for (i = pairCount - 1; i >= 0 && pairs[i].first > xMid; i--) {
            pair<Point, Point> pair = make_pair(points[pairs[i].second.first], points[pairs[i].second.second]);
            if (length(pair.first - pointMid) < length(pair.second - pointMid)) {
                isDiscarded[pairs[i].second.first] = true;
            }
            else {
                isDiscarded[pairs[i].second.second] = true;
            }
        }
        for (int i = 0; i < n; i++) if (!isDiscarded[i]) filter[s++] = points[i];
    }

    return constraintYCoordianteCenter(filter, s, y);
}

int main() {
#ifndef ONLINE_JUDGE
    freopen("inputs/input.txt", "r", stdin);
#endif

    Point p[N];
    int n;
    double y;
    scanf("%d%lf", &n, &y);
    for (int i = 0; i < n; i++) {
        scanf("%lf%lf", &p[i].x, &p[i].y);
    }
    pair<double, double> ans = constraintYCoordianteCenter(p, n, y);
    printf("x = %.5f  r = %.5f\n", ans.first, ans.second);
    return 0;
//    printf("%.5f\n", getPerpendicularBisectorIntersection(Point(1, 2), Point(2, 3), X_AXIS, 0));
}