#include <omp.h>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include "excerpt.h"

#pragma ide diagnostic ignored "openmp-use-default-none"

#define MAX_DISTANCE 10e-5

using namespace std;
enum {
    x1 = 0, x2, x3
};

typedef float fp_t;

complex<fp_t> fmaComplex(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c) {
    return {pr_product_difference(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
            fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))};
}

#define SCF(x) static_cast<fp_t>(x)
#define SCCF(x) static_cast<complex<fp_t>>(x)

template<typename fp_t>
void solve(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots) {
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

//    cout << "x^3+" << B << "x^2+" << C << "x+" << D << "=0" << endl;

    if (abs(A) < numeric_limits<fp_t>::epsilon()) {
        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
        auto p = -C / SCF(2.0L);
        auto q = sqrt(SCCF(pr_product_difference(p, p, B, D)));
        if (abs(q.imag()) < numeric_limits<fp_t>::epsilon()) {
            auto r = p + copysign(q.real(), p);
            if (abs(r) < numeric_limits<fp_t>::epsilon()) {
                roots[x1] = D / B;
                roots[x2] = -roots[x1];
            } else {
                roots[x1] = D / r;
                roots[x2] = r / B;
            }
        } else {
            roots[x1] = (p + q) / B;
            roots[x2] = (p - q) / B;
        }
    } else {
        auto b = -B / (A * SCF(3.0L));
        auto c = C / A;
        auto d = D / A;
        auto s = fma(SCF(3.0L) * b, b, -c);
        auto t = fma((fma(-b, b, s)), b, -d);
        complex<fp_t> y1, y2;
        if (abs(s) < numeric_limits<fp_t>::epsilon()) {
            if (abs(t) < numeric_limits<fp_t>::epsilon()) {
                t = copysign(numeric_limits<fp_t>::epsilon(), t) + t;
            }
//            y1 = pow(SCCF(-t), SCF(1.0L / 3.0L));
            y1 = cbrt(-t);
            y2 = y1 * (complex<fp_t>(SCF(-1.0L), sqrt(SCF(3.0L)))) / SCF(2.0L);
            roots[x1] = b - y1;
            roots[x2] = fmaComplex(-y1, complex<fp_t>(SCF(-1.0L), sqrt(SCF(3.0L))) / SCF(2.0L), b);
        } else {
            auto u = sqrt((SCF(4.0L / 3.0L) * s));
            auto v = asin(SCCF(SCF(3.0L) * t / s) / u) / SCF(3.0L);
            auto w = copysign(numbers::pi_v<fp_t> / SCF(3.0L), v.real()) - v;
            y1 = u * sin(v);
            y2 = u * sin(w);
            roots[x1] = fmaComplex(-u, sin(v), b);
            roots[x2] = fmaComplex(-u, sin(w), b);
        }
        roots[x3] = y1 + y2 + b;
    }

}

template<typename fp_t>
void comparator(vector<fp_t> &rootsTruth, vector<fp_t> &rootsOut, fp_t &absOut, fp_t &relOut) {
    double abs = 10000.0;
    double rel = 10000.0;
    for (int i = 0; i < rootsOut.size(); i++) {
        double absLoc = std::abs(double(rootsTruth[i]) - double(rootsOut[i]));
        abs = min(absLoc, abs);
        rel = min(std::abs(
                double(absLoc + std::numeric_limits<fp_t>::epsilon()) /
                double(max(rootsOut[i], rootsTruth[i]) + std::numeric_limits<fp_t>::epsilon())), rel);
    }
    absOut = abs;
    relOut = rel;
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count, vector<fp_t> &coeffs) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0,
                              MAX_DISTANCE, -1, 1, roots, coefficients);
    vector<complex<fp_t>> roots_computed(roots_count);
    solve<fp_t>(coefficients, roots_computed);
//    compare_roots_complex<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
//                                max_absolute_error, max_relative_error);

    std::vector<fp_t> roots_to_check_parsed;
    for (auto root: roots_computed) {
        if (std::numeric_limits<fp_t>::epsilon() > abs(root.imag())) {
            roots_to_check_parsed.push_back(root.real());
        }
    }

    comparator<fp_t>(roots, roots_to_check_parsed, max_absolute_error, max_relative_error);
    coeffs = coefficients;
    return pair<fp_t, fp_t>(max_absolute_error, max_relative_error);
}

int main() {
    fp_t max_deviation_abs = 0, max_deviation_rel = 0;
    auto cores = omp_get_num_procs();
    auto *deviations_abs = new fp_t[cores];
    auto *deviations_rel = new fp_t[cores];
    //test2();
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

//int main() {
//    unsigned int roots = 3;
//    vector<fp_t> roots_computed(roots);
//    auto result = testPolynomial<fp_t>(roots, roots_computed);
//    cout << "Max deviation absolute: " << result.first << endl;
//    cout << "Max deviation relative: " << result.second << endl;
//}