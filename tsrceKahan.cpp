//
// Created by eszdman on 12.12.22.
//
#include "tsrceKahan.h"
#include "excerpt.h"

#define gen(fp_t2) template void tsrceKahanReal(vector<fp_t2> &coefficients, vector<fp_t2> &roots); \
template void tsrceKahan(vector<fp_t2> &coefficients, vector<complex<fp_t2>> &roots);

template<typename fp_t> complex<fp_t> fmaComplex(complex<fp_t> a, complex<fp_t> b, complex<fp_t> c) {
    return {pr_product_difference(a.real(), b.real(), a.imag(), b.imag()) + c.real(),
            fma(a.real(), b.imag(), fma(a.imag(), b.real(), c.imag()))};
}

template<typename fp_t>
void tsrceKahanReal(vector<fp_t> &coefficients, vector<fp_t> &roots){
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

    //cout << "x^3+" << B << "x^2+" << C << "x+" << D << "=0" << endl;

    if (abs(A) < numeric_limits<fp_t>::epsilon()) {
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
            roots[x1] = (p + q) / B;
            roots[x2] = (p - q) / B;
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
            y2 = fmaComplex<fp_t>(y1, -fmaComplex<fp_t>((1.iF), sqrt(SCFC(3.0L)), SCFC(-1.0L)) / SCFC(2.0L), b);
            y1 = b - y1;
        } else {
            auto u = sqrt(SCFC(s * SCF(4.0L / 3.0L)));
            auto v = asin(SCFC(SCF(3.0L) * t / (s * u))) / SCF(3.0L);
            auto w = copysign((numbers::pi_v<fp_t> / SCF(3.0L)), v.real()) - v;
            y1 = fmaComplex<fp_t>(-sin(v), u, b);
            y2 = fmaComplex<fp_t>(-sin(w), u, b);

        }
        roots[x1] = copysign(sqrt((y1*complex<fp_t>(y1.real(),-y1.imag())).real()),y1.real());
        roots[x2] = copysign(sqrt((y2*complex<fp_t>(y2.real(),-y2.imag())).real()),y2.real());
        roots[x3] = roots[x1] + roots[x2] + b;
        /*if(isnan(roots[x3])){
            cout<<"Nan"<<endl;
        }*/
    }
}
template<typename fp_t>
void tsrceKahan(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots){
    auto A = coefficients[3];
    auto B = coefficients[2];
    auto C = coefficients[1];
    auto D = coefficients[0];

//    cout << "x^3+" << B << "x^2+" << C << "x+" << D << "=0" << endl;

    if (abs(A) < numeric_limits<fp_t>::epsilon()) {
        roots[x3] = (abs(B) + abs(C) + abs(D)) / A;
        auto p = -C / SCF(2.0L);
        auto q = sqrt(SCFC(pr_product_difference(p, p, B, D)));
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
            y1 = cbrt(-t);
            roots[x1] = b - y1;
            roots[x2] = fmaComplex<fp_t>(y1, -SCFC((SCFC(1.iF) * SCFC(sqrt(SCF(3.0L))) - SCFC(1.0L)) / SCFC(2.0L)), b);
            roots[x3] = y1 + fmaComplex<fp_t>(y1, SCFC((SCFC(1.iF) * SCFC(sqrt(SCF(3.0L))) - SCFC(1.0L)) / SCFC(2.0L)), b);
            return;
        } else {
            auto u = sqrt(SCFC(s / SCF(4.0L / 3.0L)));
            auto v = asin(SCFC(SCF(3.0L) * t / (s * u))) / SCF(3.0L);
            auto w = copysign((numbers::pi_v<fp_t> / SCF(3.0L)), v.real()) - v;
            roots[x1] = fmaComplex<fp_t>(-sin(v), u, b);
            roots[x2] = fmaComplex<fp_t>(-sin(w), u, b);
            roots[x3] = roots[x1] + roots[x2] + SCFC(b);
            return;
        }
    }
}

gen(float)
gen(double)
gen(long double)