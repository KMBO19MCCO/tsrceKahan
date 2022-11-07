#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "excerpt.h"

using namespace std;
enum {
    x1 = 0, x2, x3
};
typedef float fp_t;


template<typename fp_t>
int solve(vector<fp_t> &coefficients, vector<fp_t> &roots) {
    fp_t A = coefficients[3];
    fp_t B = coefficients[2];
    fp_t C = coefficients[1];
    fp_t D = coefficients[0];
    roots[x3] = A != 0 ? (abs(B) + abs(C) + abs(D)) / A : 0;
    auto p = C / static_cast<fp_t>(-2L);
    auto q = sqrt(fma(p, p, -B * D));
    if (q >= 0) {
        auto r = p + copysign(1.f, p) * q;
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

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t deviation;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, numeric_limits<fp_t>::min(), -1, 1, roots, coefficients);
    vector<fp_t> roots_computed(roots_count);
    solve<fp_t>(coefficients, roots_computed);
    auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, deviation);
    switch (result) {
        case PR_2_INFINITE_ROOTS:
            cout << "INFINITE ROOTS";
            break;
        case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
            cout << "AT LEAST ONE ROOT IS FAKE";
            break;
        case PR_AT_LEAST_ONE_ROOT_LOST:
            cout << "AT LEAST ONE ROOT LOST";
            break;
        default:
            break;
    }
    return deviation;
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
