#define __FMA__ 1

#include <omp.h>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <excerpt.h>
#include "tsrceKahan.h"

#pragma ide diagnostic ignored "openmp-use-default-none"

#define MAX_DISTANCE 1e-25

using namespace std;

typedef double fp_t;



template<typename fp_t>
auto testPolynomial(unsigned int roots_count, vector<fp_t> &coeffs) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0,
                              MAX_DISTANCE, MAX_DISTANCE, 1e+1, roots, coefficients);
    vector<fp_t> roots_computed(roots_count);
    if (roots_count + 1 < 4) {
        coefficients.resize(4);
        coefficients[3] = 0;
    }
    tsrceKahanReal<fp_t>(coefficients, roots_computed);
    //compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
    //                    max_absolute_error, max_relative_error);
    compare_roots2<fp_t>(roots.size(),roots.size(),roots, roots_computed, max_absolute_error, max_relative_error);
    coeffs = coefficients;
    return pair<fp_t, fp_t>(max_absolute_error, max_relative_error);
}

void test2() {
    unsigned int roots_count = 3;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    vector<complex<fp_t>> roots_computed(roots_count);
    coefficients[0] = -0.999957;
    coefficients[1] = 2.99991;
    coefficients[2] = -2.99996;
    coefficients[3] = 1;
    tsrceKahan<fp_t>(coefficients, roots_computed);
    cout << roots_computed[0] << " " << roots_computed[1] << " " << roots_computed[2];
    exit(0);
}


int main() {
    //test2();
    fp_t max_deviation_abs = 0, max_deviation_rel = 0;
    auto cores = omp_get_num_procs();
    auto *deviations_abs = new fp_t[cores];
    auto *deviations_rel = new fp_t[cores];
#pragma omp parallel for
    for (auto i = 0; i < cores; ++i) {
        deviations_abs[i] = 0;
        deviations_rel[i] = 0;
    }

    cout << "Threads: " << cores << endl;
    cout << "Started" << endl;
//    auto timer = Timer("compute");
//    timer.start();
    int roots = 3;

#pragma omp parallel for
    for (auto i = 0; i < 1000'1000; ++i) {
        auto thread_id = omp_get_thread_num();
        vector<fp_t> roots_computed(roots);
        auto result = testPolynomial<fp_t>(roots, roots_computed);
        if (result.first > deviations_abs[thread_id] and !isinf(result.first)) {
            deviations_abs[thread_id] = result.first;
        }
        if (result.second > deviations_rel[thread_id] and !isinf(result.second)) {
            deviations_rel[thread_id] = result.second;
        }
    }

    cout << "Computed, started searching" << endl;
//    timer.stop();
    for (auto i = 0; i < cores; ++i) {
        if (deviations_abs[i] > max_deviation_abs) {
            max_deviation_abs = deviations_abs[i];
        }
        if (deviations_rel[i] > max_deviation_rel) {
            max_deviation_rel = deviations_rel[i];
        }
    }
    cout << "Max deviation absolute: " << max_deviation_abs << endl;
    cout << "Max deviation relative: " << max_deviation_rel << endl;
    delete[] deviations_abs;
    delete[] deviations_rel;
    return 0;
}
