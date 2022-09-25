#include <iostream>
#include <vector>
#include <cmath>
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

    solution.at(x3) =
            coefficients.at(A) != 0 ? (abs(coefficients.at(B)) + abs(coefficients.at(C)) + abs(coefficients.at(D))) /
                                      coefficients.at(A) : 0;
    auto p = -coefficients.at(C) / 2;
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

    return -1;
}

int main() {
    vector<float> coefficients = {0, 2, -6, 1};
    vector<float> solution = {.0, .0, .0};
    cout << solve(coefficients, solution) << endl;
    for (auto x: solution) {
        cout << x << ' ';
    }
    return 0;
}
