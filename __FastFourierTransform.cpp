//
// Created by Armeria on 2020/3/13.
//

#include <cstdio>
#include <cmath>
#include <algorithm>
// include <complex>

using namespace std;

//typedef complex<double> Complex;

struct Complex {
    double r;
    double i;

    Complex() = default;

    Complex(double r, double i) {
        this->r = r;
        this->i = i;
    }

    double real() { return r; }
    double imaginary() { return i; }

    Complex operator+(const Complex &A) const {
        Complex ret;
        ret.r = (r + A.r);
        ret.i = (i + A.i);
        return ret;
    }

    Complex operator-(const Complex &A) const {
        Complex ret;
        ret.r = (r - A.r);
        ret.i = (i - A.i);
        return ret;
    }

    Complex operator*(const Complex &A) const {
        Complex ret;
        ret.r = (r * A.r - i * A.i);
        ret.i = (r * A.i + i * A.r);
        return ret;
    }
};

const int N = (1 << 21) + 5;
const Complex I = Complex(0, 1);
const double PI = acos(-1.0);

Complex a[N], b[N], res[N << 1];
double x;
int n, m, s;

void butterfly(Complex f[], int n) {
    for (int i = 1, j = (n >> 1); i < n - 1; i++) {
        if (i < j) swap(f[i], f[j]);
        int k = n >> 1;
        while (j >= k) {
            j = j - k;
            k >>= 1;
        }
        if (j < k) j += k;
    }
}

void fft(Complex f[], int n, int rev) {
    butterfly(f, n);
    for (int len = 2; len <= n; len <<= 1) {
        Complex step = Complex(cos(2 * PI / len), sin(rev * 2 * PI / len));
//        Complex step = exp(I * (2.0 * M_PI / length * rev));
        for (int j = 0; j < n; j += len) {
            Complex cur = Complex(1, 0);
            for (int k = j; k < j + (len >> 1); k++) {
                Complex x = f[k];
                Complex y = cur * f[k + (len >> 1)];
                f[k] = x + y;
                f[k + (len >> 1)] = x - y;
                cur = cur * step;
            }
        }
    }
}

Complex tmp[N];

void dft(Complex *f, int n, int rev) {
    if (n == 1) return;
    for (int i = 0; i < n; i++) tmp[i] = f[i];
    for (int i = 0; i < n; i++)
        if (i & 1)
            f[(n >> 1) + (i >> 1)] = tmp[i];
        else
            f[(i >> 1)] = tmp[i];

    Complex *g = f, *h = f + (n >> 1);
    dft(g, n >> 1, rev);
    dft(h, n >> 1, rev);

    Complex cur = Complex(1, 0);
    Complex step = Complex(cos(2 * PI / n), sin(rev * 2 * PI / n));
    for (int i = 0; i < (n >> 1); i++) {
        tmp[i] = g[i] + cur * h[i];
        tmp[i + (n >> 1)] = g[i] - cur * h[i];
        cur = cur * step;
    }
    for (int i = 0; i < n; i++) f[i] = tmp[i];
}

void expand(Complex f[], int &n, int target) {
    int pow2 = 1;
    while (target > pow2) pow2 <<= 1;
    for (int i = n; i < pow2; i++) f[i] = Complex(0, 0);
    n = pow2;
}

void solveRecursive() {
    dft(a, n, 1);
    dft(b, m, 1);

    for (int i = 0; i < n; i++) {
        res[i] = a[i] * b[i];
    }

    dft(res, n, -1);
}

void solveNonRecursive() {
    fft(a, n, 1);
    fft(b, m, 1);

    for (int i = 0; i < n; i++) {
        res[i] = a[i] * b[i];
    }

    fft(res, n, -1);
}

int main() {

#ifndef ONLINE_JUDGE
    freopen("inputs/input.txt", "r", stdin);
#endif

    scanf("%d%d", &n, &m);
    n++, m++;
    s = n + m - 1;
    for (int i = 0; i < n; i++) scanf("%lf", &x), a[i] = Complex(x, 0);
    for (int i = 0; i < m; i++) scanf("%lf", &x), b[i] = Complex(x, 0);

    int target = ((max(n, m)) << 1);
    expand(a, n, target);
    expand(b, m, target);

//    solveRecursive();
    solveNonRecursive();

    for (int i = 0; i < s; i++) {
        printf("%d ", (int) (res[i].real() / n + 0.5));
    }

    return 0;
}
