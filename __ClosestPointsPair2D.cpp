//
// Created by Armeria on 2020/3/14.
//

#include <cstdio>
#include <algorithm>
#include <cmath>

using namespace std;

const int N = 200000 + 5, C = 10;

struct Point {
    int x, y;
} points[N];

int n;
double ans;
pair<Point, Point> cp;

double cmp_x(Point a, Point b) {
    return a.x == b.x ? a.y < b.y : a.x < b.x;
}

double cmp_y(Point a, Point b) {
    return a.y < b.y;
}

void updateAnswer(const Point a, const Point b) {
    double dist = sqrt(1.0 * (a.x - b.x) * (a.x - b.x) + 1.0 * (a.y - b.y) * (a.y - b.y));
//    printf("%.5lf\n", dist);
    if (dist < ans) {
        ans = dist;
        cp = make_pair(a, b);
    }
}

void findClosest(Point points[], int n) {
    if (n <= 3) {
        for (int i = 0; i < n - 1; i++)
            for (int j = i + 1; j < n; j++)
                updateAnswer(points[i], points[j]);
        sort(points, points + n, cmp_y);
        return;
    }

    int k = (n >> 1);
    int xMid = points[k].x;
    findClosest(points, k);
    findClosest(points + k, n - k);

    static Point tmp[N];
    merge(points, points + k, points + k, points + n, tmp, cmp_y);
    copy(tmp, tmp + n, points);

    int cur = 0;
    static Point legal[C];
    for (int i = 0; i < n; i++) {
        if (abs(points[i].x - xMid) < ans) {
            for (int j = cur - 1; j >= 0; j--) {
                if (points[i].y - legal[j].y > ans) break;
                updateAnswer(points[i], legal[j]);
            }
            legal[cur++] = points[i];
        }
    }

}

int main() {

#ifndef ONLINE_JUDGE
//    freopen("inputs/input.txt", "r", stdin);
    freopen("inputs/P1429_1.in", "r", stdin);
#endif

    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        scanf("%d%d", &points[i].x, &points[i].y);
    }
    sort(points, points + n, cmp_x);

    ans = 1e17;
    findClosest(points, n);
    printf("%.4f\n", ans);

    return 0;
}
