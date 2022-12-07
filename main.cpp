#include <omp.h>
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <excerpt.h>

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
#define SCFC(x) static_cast<complex<fp_t>>(x)

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
        auto p = -C / SCF(2.0L);
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
        auto b = -B / (A * SCF(3.0L));
        auto c = C / A;
        auto d = D / A;
        auto s = fma(SCF(3.0L) * b, b, -c);
        //auto s = pr_product_difference(static_cast<fp_t>(3.0L) * b, b, c, static_cast<fp_t>(1.0L));
        auto t = fma(fma(-b, b, s), b, -d);
        //auto t = pr_product_difference(pr_product_difference(static_cast<fp_t>(1.0L), s, b, b), b, d,
        //                               static_cast<fp_t>(1.0L));
        //auto t = pr_product_difference(b, s, b, b * b) - d;
        complex<fp_t> y1, y2;
        if (abs(s) < numeric_limits<fp_t>::epsilon()) {
            y1 = cbrt(-t);
            roots[x1] = b - y1;
            roots[x2] = fmaComplex(y1, -SCFC((SCFC(1.iF) * SCFC(sqrt(SCF(3.0L))) - SCFC(1.0L)) / SCFC(2.0L)), b);
            roots[x3] = y1 + fmaComplex(y1, SCFC((SCFC(1.iF) * SCFC(sqrt(SCF(3.0L))) - SCFC(1.0L)) / SCFC(2.0L)), b);
            return;
        } else {
            auto u = sqrt(SCFC(s / SCF(4.0L / 3.0L)));
            auto v = asin(SCFC(SCF(3.0L) * t / (s * u))) / SCF(3.0L);
            auto w = copysign((numbers::pi_v<fp_t> / SCF(3.0L)), v.real()) - v;
            roots[x1] = fmaComplex(-sin(v), u, b);
            roots[x2] = fmaComplex(-sin(w), u, b);
            roots[x3] = roots[x1] + roots[x2] + SCFC(b);
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

    //cout << "x^3+" << B << "x^2+" << C << "x+" << D << "=0" << endl;

    if (abs(A) < numeric_limits<fp_t>::epsilon()) { // unused
        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
        auto p = -C / SCF(2.0L);
        auto q = sqrt(pr_product_difference(p, p, B, D));
        if (abs(q) < numeric_limits<fp_t>::epsilon()) {
            auto r = p + copysign(q, p * q);
            if (abs(r) < numeric_limits<fp_t>::epsilon()) {
                roots[x1] = D / B;
                roots[x2] = -roots[x1];
            } else {
                roots[x1] = D / r;
                roots[x2] = r / B;
            }
        } else {
            roots[x1] = p + q / B;
            roots[x2] = p - q / B;
        }
    } else {
        auto b = -B / (A * SCF(3.0L));
        auto c = C / A;
        auto d = D / A;
        auto s = fma(SCF(3.0L) * b, b, -c);
        auto t = fma(fma(-b, b, s), b, -d);
        complex<fp_t> y1, y2;
        if (abs(s) < numeric_limits<fp_t>::epsilon()) {
            y1 = cbrt(-t);
            y2 = fmaComplex(y1, -SCFC((SCFC(1.iF) * SCFC(sqrt(SCF(3.0L))) - SCFC(1.0L)) / SCFC(2.0L)), b);
            y1 = b - y1;
        } else {
            auto u = sqrt(s * SCF(4.0L / 3.0L));
            auto v = asin(SCFC(SCF(3.0L) * t / (s * u))) / SCF(3.0L);
            auto w = copysign((numbers::pi_v<fp_t> / SCF(3.0L)), v.real()) - v;
            y1 = fmaComplex(-sin(v), u, b);
            y2 = fmaComplex(-sin(w), u, b);
        }
        //roots[x1] = copysign(hypot(y1.real(), y1.imag()), y1.real());
        //roots[x2] = copysign(hypot(y2.real(), y2.imag()), y2.real());
        //roots[x1] = copysign(sqrt((y1*complex<fp_t>(y1.real(),-y1.imag())).real()),y1.real());
        //roots[x2] = copysign(sqrt((y2*complex<fp_t>(y2.real(),-y2.imag())).real()),y2.real());
        //roots[x1] = (y1*complex<fp_t>(0,-y1.imag())).real();
        //roots[x2] = (y2*complex<fp_t>(0,-y2.imag())).real();
        roots[x1] = y1.real();
        roots[x2] = y2.real();
        roots[x3] = roots[x1] + roots[x2] + b;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count, vector<fp_t> &coeffs) {
    fp_t max_absolute_error, max_relative_error;
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0,
                              MAX_DISTANCE, -1, 1, roots, coefficients);
    vector<fp_t> roots_computed(roots_count);
    if (roots_count + 1 < 4) {
        coefficients.resize(4);
        coefficients[3] = 0;
    }
    solveReal<fp_t>(coefficients, roots_computed);
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
    solve<fp_t>(coefficients, roots_computed);
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
