#include <omp.h>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include "excerpt.h"

#define MAX_DISTANCE 10e-5

using namespace std;
enum {
    x1 = 0, x2, x3
};
typedef float fp_t;


template<typename fp_t>
void solve(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots) {
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

    if (std::numeric_limits<fp_t>::epsilon() > abs(A)) {
//        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
//        roots[x3] = std::numeric_limits<fp_t>::epsilon(); // maybe framework bug
        roots[x3] = 0; // maybe framework bug
        auto p = -C / 2;
        auto q = sqrt(fma(p, p, -B * D));
        if (std::numeric_limits<fp_t>::epsilon() > abs(q)) {
            auto r = fma(copysign(1, q), q, p);
            if (std::numeric_limits<fp_t>::epsilon() > abs(r)) {
                roots[x1] = D / B;
                roots[x2] = -roots[x1];
                return;
            } else {
                roots[x1] = D / r;
                roots[x2] = r / B;
                return;
            }
        } else {
            roots[x1] = p + q / B;
//            roots[x1] = fma(static_cast<fp_t>(1LL) / B, p, p);
            roots[x2] = p - q / B;
//            roots[x2] = fma(static_cast<fp_t>(1LL) / B, -p, p);
            return;
        }
    } else {
        auto b = -B / (A * 3);
        auto c = C / A;
        auto d = D / A;
        auto s = 3 * fma(b, b, -c);
        auto t = fma(fma(-b, b, s), b, -d);
        complex<fp_t> y1, y2;
        if (s == 0) {
            y1 = pow(-t, static_cast<fp_t>(1.0) / 3.0);
            y2 = y1 * static_cast<fp_t>((1.iF * sqrt(3) - 1) / 2);
        } else {
            auto u = sqrt(static_cast<complex<fp_t>>(s / 3) * static_cast<fp_t>(4));
            auto v = asin(static_cast<complex<fp_t>>(t) / (s * u));
            auto w = (numbers::pi_v<fp_t> / 3) - v;
            y1 = sin(v) * u;
            y2 = sin(w) * u;
        }
        roots[x1] = b - y1;
        roots[x2] = b - y2;
        roots[x3] = y1 + y2 + b;
        return;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count, vector<fp_t> &coeffs) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots, coefficients);
    vector<complex<fp_t>> roots_computed(roots_count);
    solve<fp_t>(coefficients, roots_computed);
    compare_roots_complex<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                                max_absolute_error, max_relative_error);
    coeffs = coefficients;
    return max_absolute_error;
}

int main() {
    fp_t max_deviation = 0;
    auto cores = omp_get_num_procs();
    auto *deviations = new fp_t[cores];

#pragma omp parallel for
    for (auto i = 0; i < cores; ++i) {
        deviations[i] = 0;
    }

    cout << "Threads: " << cores << endl;
    cout << "Started" << endl;
//    auto timer = Timer("compute");
//    timer.start();

#pragma omp parallel for
    for (auto i = 0; i < 1000'1000; ++i) {
        auto thread_id = omp_get_thread_num();
        vector<fp_t> roots_computed(4);
        auto deviation = testPolynomial<fp_t>(3, roots_computed);
        if (deviation > deviations[thread_id]) {
            deviations[thread_id] = deviation;
            for (auto root: roots_computed){
                cout << root << ' ';
            }
            cout << endl;
        }
    }

    cout << "Computed, started searching" << endl;
//    timer.stop();
    for (auto i = 0; i < cores; ++i) {
        if (deviations[i] > max_deviation) {
            max_deviation = deviations[i];
        }
    }
    cout << "Max deviation: " << max_deviation << endl;
    delete[] deviations;
    return 0;
}
