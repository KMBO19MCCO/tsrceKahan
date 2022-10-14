#include <iostream>
#include <vector>
#include <cmath>
#include "excerpt.h"
using namespace std;
enum {
    x1 = 0, x2, x3
};
typedef float fp_t;
template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template<typename fp_t>
int solve(vector<fp_t> &coeff, vector<fp_t> &roots) {
    fp_t A = coeff[3];
    fp_t B = coeff[2];
    fp_t C = coeff[1];
    fp_t D = coeff[0];
    roots[x3] = A != 0 ? (abs(B) + abs(C) + abs(D))/A : 0;
    auto p = C / static_cast<fp_t>(-2L);
    auto q = sqrt(p * p - B * D);
    if (q >= 0) {
        auto r = p + sgn(p) * q;
        if (r == 0) {
            roots[x1] = D / B;
            roots[x2] = -x1;
            return 0;
        } else {
            roots[x1] = D / r;
            roots[x2] = r / B;
            return 1;
        }
    } else {
        roots[x1] = p / B + q / B;
        roots[x2] = p / B - q / B;
        return 2;
    }
}
int main() {
    int p = 3;
    vector<fp_t> coefficients(p+1);
    vector<fp_t> roots(p);
    vector<fp_t> output(p);
    generate_polynomial<fp_t>(p, 0, p, 0, 8.0L/5, -5, 5.0L, roots, coefficients);
    fp_t dev;
    cout << solve(coefficients, output) << endl;
    cout<<compare_roots<fp_t>(p,p,output,roots,dev)<<endl;
    cout<<"dev:"<<dev<<endl<<"roots:";
    for (float root : roots) {
        cout<<setprecision(numeric_limits<fp_t>::digits10 + 1)<< root << ' ';
    }
    cout<<endl<<"output:";
    for (float root : output) {
        cout<<setprecision(numeric_limits<fp_t>::digits10 + 1)<< root << ' ';
    }

    return 0;
}
