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
complex<fp_t> fmaComplex(complex<fp_t> a,complex<fp_t> b,complex<fp_t> c){
    return {pr_product_difference(a.real(),b.real(),a.imag(),b.imag()) + c.real(),
                         fma(a.real(),b.imag(),fma(a.imag(),b.real(),c.imag()))};
}


template<typename fp_t>
void solve(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots) {
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

    if (std::numeric_limits<fp_t>::epsilon() > abs(A)) { // unused
//        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
//        roots[x3] = std::numeric_limits<fp_t>::epsilon(); // maybe framework bug
        roots[x3] = 0; // maybe framework bug
        auto p = -C / static_cast<fp_t>(2.0L);
        auto q = sqrt(pr_product_difference(p, p, B, D));
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
            //roots[x1] = fma(static_cast<fp_t>(1.0L) / B, q, p);
            roots[x2] = p - q / B;
            //roots[x2] = fma(static_cast<fp_t>(1.0L) / B, -q, p);
            return;
        }
    } else {
        auto b = -B / (A * static_cast<fp_t>(3.0L));
        auto c = C / A;
        auto d = D / A;
//        auto s = fma(static_cast<fp_t>(3.0L) * b, b, -c);
        auto s = pr_product_difference(static_cast<fp_t>(3.0L) * b, b, c, static_cast<fp_t>(1.0L));
//        auto t = fma(fma(-b, b, s), b, -d);
        auto t = pr_product_difference(pr_product_difference(static_cast<fp_t>(1.0L), s, b, b), b, d,
                                       static_cast<fp_t>(1.0L));
//        auto t = pr_product_difference(b, s, b, b * b) - d;
        complex<fp_t> y1, y2;
        if (s == 0) {
            y1 = pow(-t, static_cast<fp_t>(1.0L / 3.0L));
            y2 = y1 * static_cast<complex<fp_t>>((static_cast<complex<fp_t>>(1.iF) *
                                                  static_cast<complex<fp_t>>(sqrt(static_cast<fp_t>(3.0L))) -
                                                  static_cast<complex<fp_t>>(1.0L))
                                                 / static_cast<complex<fp_t>>(2.0L));
            roots[x1] = b - y1;
            roots[x2] = b - y2;
            roots[x3] = y1 + y2 + b;
            return;
        } else {
            auto u = sqrt(static_cast<complex<fp_t>>(s / static_cast<fp_t>(3.0L) * static_cast<fp_t>(4.0L)));
            auto v = asin(static_cast<complex<fp_t>>(static_cast<fp_t>(3.0L) * t) / (s * u)) / static_cast<fp_t>(3.0L);
            auto w = (numbers::pi_v<fp_t> / static_cast<fp_t>(3.0L)) * copysign(static_cast<fp_t>(1.0L), v.real()) - v;
            roots[x1] = fmaComplex(-sin(v),u,b);
            roots[x2] = fmaComplex(-sin(w),u,b);
            roots[x3] = roots[x1] + roots[x2]  + b;
            return;
        }
    }
}
template<typename fp_t>
void solveReal(vector<fp_t> &coefficients, vector<fp_t> &roots) {
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

    if (std::numeric_limits<fp_t>::epsilon() > abs(A)) { // unused
//        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
//        roots[x3] = std::numeric_limits<fp_t>::epsilon(); // maybe framework bug
        roots[x3] = 0; // maybe framework bug
        auto p = -C / static_cast<fp_t>(2.0L);
        auto q = sqrt(pr_product_difference(p, p, B, D));
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
            //roots[x1] = fma(static_cast<fp_t>(1.0L) / B, q, p);
            roots[x2] = p - q / B;
            //roots[x2] = fma(static_cast<fp_t>(1.0L) / B, -q, p);
            return;
        }
    } else {
        auto b = -B / (A * static_cast<fp_t>(3.0L));
        auto c = C / A;
        auto d = D / A;
//        auto s = fma(static_cast<fp_t>(3.0L) * b, b, -c);
        auto s = pr_product_difference(static_cast<fp_t>(3.0L) * b, b, c, static_cast<fp_t>(1.0L));
//        auto t = fma(fma(-b, b, s), b, -d);
        auto t = pr_product_difference(pr_product_difference(static_cast<fp_t>(1.0L), s, b, b), b, d,
                                       static_cast<fp_t>(1.0L));
//        auto t = pr_product_difference(b, s, b, b * b) - d;
        complex<fp_t> y1, y2;
        if (s == 0) {
            y1 = pow(-t, static_cast<fp_t>(1.0L / 3.0L));
            y2 = y1 * static_cast<complex<fp_t>>((static_cast<complex<fp_t>>(1.iF) *
                                                  static_cast<complex<fp_t>>(
                                                          sqrt(static_cast<fp_t>(3.0L))) -
                                                  static_cast<complex<fp_t>>(1.0L))
                                                 / static_cast<complex<fp_t>>(2.0L));
            y1 = b - y1;
            y2 = b - y2;
        } else {
            auto u = sqrt(static_cast<complex<fp_t>>(s / static_cast<fp_t>(3.0L) * static_cast<fp_t>(4.0L)));
            auto v = asin(static_cast<complex<fp_t>>(static_cast<fp_t>(3.0L) * t) / (s * u)) / static_cast<fp_t>(3.0L);
            auto w = (numbers::pi_v<fp_t> / static_cast<fp_t>(3.0L)) * copysign(static_cast<fp_t>(1.0L), v.real()) - v;
            y1 = fmaComplex(-sin(v),u,b);
            y2 = fmaComplex(-sin(w),u,b);
        }
        roots[x1] = copysign(hypot(y1.real(),y1.imag()),y1.real());
        roots[x2] = copysign(hypot(y2.real(),y2.imag()),y2.real());
        roots[x3] = roots[x1] + roots[x2]  + b;
        return;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count, vector<fp_t> &coeffs) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, MAX_DISTANCE, -1, 1, roots, coefficients);
    vector<fp_t> roots_computed(roots_count);
    solveReal<fp_t>(coefficients, roots_computed);
    compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots,
                                max_absolute_error, max_relative_error);
    coeffs = coefficients;
    return max_absolute_error;
}

void test2() {
    unsigned int roots_count = 3;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    vector<complex<fp_t>> roots_computed(roots_count);
    coefficients[0] = -0.999957;
    coefficients[1] = 2.99991;
    coefficients[2] = -2.99996;
    coefficients[3] = 1;
    solve<fp_t>(coefficients, roots_computed);
    exit(-1);
}

int main() {
    fp_t max_deviation = 0;
    auto cores = omp_get_num_procs();
    auto *deviations = new fp_t[cores];
    //test2();
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
        if (deviation > deviations[thread_id] and !isinf(deviation)) {
            deviations[thread_id] = deviation;
            //cout<<deviation<<endl;
//            for (auto root: roots_computed){
//                cout << root << ' ';
//            }
//            cout << endl;
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
