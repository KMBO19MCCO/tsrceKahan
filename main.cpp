#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include "excerpt.h"

using namespace std;

enum {
    A = 0, B, C, D
};

enum {
    x1 = 0, x2, x3
};

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename fp_t>
int solve(vector<fp_t> &coefficients, vector<fp_t> &solution) {
    if (coefficients.at(A) == 0) {
        solution.at(x3) =
                (abs(coefficients.at(B)) + abs(coefficients.at(C)) + abs(coefficients.at(D))) / coefficients.at(A);
        auto p = -coefficients.at(C) / 2.0;
        auto q = sqrt(p * p - coefficients.at(B) * coefficients.at(D));
        if (q >= 0) {
            auto r = p + sgn(p) * q;
            if (r == 0) {
                solution.at(x1) = coefficients.at(D) / coefficients.at(B);
                solution.at(x2) = -x1;
                return 0;
            } else {
                solution.at(x1) = coefficients.at(D) / r;
                solution.at(x2) = r / coefficients.at(B);
                return 1;
            }
        } else {
            solution.at(x1) = p / coefficients.at(B) + q / coefficients.at(B);
            solution.at(x2) = p / coefficients.at(B) - q / coefficients.at(B);
            return 2;
        }
    } else {
        auto b = -(coefficients.at(B) / coefficients.at(A)) / 3;
        auto c = coefficients.at(C) / coefficients.at(A);
        auto d = coefficients.at(D) / coefficients.at(A);
        auto s = 3 * b * b - c;
        auto t = (s - b * b) * b - d;
        fp_t y1, y2;
        if (s == 0) {
            y1 = pow(-t, 1 / 3);
            y2 = y1 * (-1 + 1i * sqrt(3)) / 2;
        } else {
            auto u = sqrt(4 * s / 3);
            auto v = asin((3 * t / s) / u) / 3;
            auto w = (numbers::pi / 3) * sgn(real(v)) - v;
            y1 = u * sin(v);
            y2 = u * sin(w);
        }
        solution.at(x1) = b - y1;
        solution.at(x2) = b - y2;
        solution.at(x3) = y1 + y2 + b;
        return 3;
    }
    return -1;
}

int main() {
    unsigned p = 3;
    vector<double> roots(p);
    vector<double> coefficients(p + 1);
    auto result = generate_polynomial<double>(p, 0, 2, 0, 10.0 / 5, -10, 10, roots, coefficients);
//    vector<complex<double>> coefficients = {2., -6, 0., 1.};
    vector<double> solution = {.0, .0, .0};
    cout << solve(coefficients, solution) << endl;
    for (auto x: roots) {
        cout << x << ' ';
    }
    cout << endl;
    for (auto x: solution) {
        cout << x << ' ';
    }
    return 0;
}
