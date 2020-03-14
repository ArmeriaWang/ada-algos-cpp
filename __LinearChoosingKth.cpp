//
// Created by Armeria on 2020/3/14.
//

#include <cstdio>
#include <algorithm>
//#include <cmath>
#include <vector>

using namespace std;

const int N = 100;

int getKthSplit(int a[], int n, int k) {
    static const int c = 5;

    if (n < c) {
        sort(a, a + n);
        return a[k];
    }

    int s = 0;
    for (int i = 0; i < n - c + 1; i += c) {
        sort(a + i, a + i + c);
        swap(a[s++], a[i + (c >> 1)]);
    }
    int xMid = getKthSplit(a, s, ((s - 1) >> 1));
    int b[N];
    int m = 0, midNum = 0;
    for (int p = 0; p < n; p++) if (a[p] != xMid) b[m++] = a[p]; else midNum++;
    int i = 0, j = m - 1;
    while (i <= j) {
        while (b[i] < xMid) i++;
        while (b[j] > xMid) j--;
        if (i < j) {
            swap(b[i], b[j]);
            i++;
            j--;
        }
    }

    int leftNum = i, rightNum = m - i;
//    int b[N];
//    for (int p = 0; p < n; p++) if (a[p] == xMid) leftNum--;
//    for (int p = 0; p < n; p++) printf("%d ", a[p]);
//    printf("\n");
//    for (int p = 0; p < m; p++) printf("%d ", b[p]);
//    printf("\n");
//    printf("n = %d  m = %d  k = %d  lN = %d  rN = %d  xMid = %d  midNum = %d\n", n, m, k, leftNum, rightNum, xMid, midNum);
//
    if (leftNum > k) return getKthSplit(b, i, k);
    else if (leftNum + midNum > k) return xMid;
    else return getKthSplit(b + leftNum, rightNum, k - midNum - leftNum);

}

int n;
int a[N], b[N];
vector<int> sorted;
vector<int> ans;

int testAll() {
    ans.clear();
    for (int i = 0; i < n; i++) {
        copy(a, a + n, b);
        ans.push_back(getKthSplit(b, n, i));
    }

    for (int i = 0; i < n; i++) {
        if (ans[i] != sorted[i]) {
            printf("WRONG ANSWER on i = %d\n", i);
            printf("std = %d  my = %d\n", sorted[i], ans[i]);
            return 1;
        }
    }
    printf("ACCEPTED\n");
    return 0;
}

int main() {

#ifndef ONLINE_JUDGE
    freopen("inputs/input_lck.txt", "r", stdin);
#endif

    scanf("%d", &n);
    sorted.clear();
    for (int i = 0; i < n; i++) {
        scanf("%d", &a[i]);
        sorted.push_back(a[i]);
    }
    sort(sorted.begin(), sorted.end());

    copy(a, a + n, b);
    printf("%d\n", getKthSplit(b, n, 10));

    testAll();

    return 0;
}
