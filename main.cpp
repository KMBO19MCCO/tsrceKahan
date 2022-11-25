#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <complex.h>
#include "excerpt.h"

#define MAX_DISTANCE 10e-5

using namespace std;
enum {
    x1 = 0, x2, x3
};
enum {
    d = 0, c, b, a
};
typedef float fp_t;


template<typename fp_t>
int solve(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots) {
    complex<fp_t> A = coefficients[3];
    complex<fp_t> B = coefficients[2];
    complex<fp_t> C = coefficients[1];
    complex<fp_t> D = coefficients[0];

    if (coefficients[a] == 0) {
        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
        auto p = -C.real() / 2;
        auto q = sqrt(p * p - B * D);
        if (std::numeric_limits<fp_t>::epsilon() > abs(q.imag())) {
            auto r = p + copysign(q.real(), 1) * q.real();
            if (r == 0) {
                roots[x1] = D / B;
                roots[x2] = -roots[x1];
                return 0;
            } else {
                roots[x1] = D.real() / r;
                roots[x2] = r / B.real();
                return 0;
            }
        } else {
            roots[x1] = p / B.real() + q / B;
            roots[x2] = p / B.real() - q / B;
            return 0;
        }
    } else {
        auto b = -(B / A).real() / 3;
        auto c = (C / A).real();
        auto d = D / A;
        auto s = 3 * b * b - c;
        auto t = (s - b * b) * b - d;
        complex<fp_t> y1, y2;
        if (s == 0) {
            y1 = pow(-t, 1.0 / 3.0);
            y2 = static_cast<complex<fp_t>>(y1) * (static_cast<complex<fp_t>>(-1) +
                                                   static_cast<complex<fp_t>>(I) *
                                                   static_cast<complex<fp_t>>(sqrt(3))) /
                 static_cast<complex<fp_t>>(2);
        } else {
            auto u = sqrt(static_cast<complex<fp_t>>(4) * s / static_cast<complex<fp_t>>(3));
            auto v = asin(static_cast<complex<fp_t>>(3) * t / s / u) / static_cast<complex<fp_t>>(3);
            auto w = static_cast<complex<fp_t>>(numbers::pi / 3) - v;
            y1 = u * sin(v);
            y2 = u * sin(w);
        }
        roots[x1] = b - y1;
        roots[x2] = b - y2;
        roots[x3] = y1 + y2 + b;
        return 0;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots, coefficients);
    vector<complex<fp_t>> roots_computed(roots_count);
    solve<fp_t>(coefficients, roots_computed);
    compare_roots_complex<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                                              max_absolute_error, max_relative_error);
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
        auto deviation = testPolynomial<fp_t>(3);
        if (deviation > deviations[thread_id] and deviation != numeric_limits<fp_t>::infinity()) {
            deviations[thread_id] = deviation;
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
