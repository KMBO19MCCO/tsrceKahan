//
// Created by eszdman on 12.12.22.
//

#ifndef TSRCEKAHAN_TSRCEKAHAN_H
#define TSRCEKAHAN_TSRCEKAHAN_H

#define SCF(x) static_cast<fp_t>(x)
#define SCFC(x) static_cast<complex<fp_t>>(x)

#include <vector>
#include <complex>

using namespace std;
enum {
    x1 = 0, x2, x3
};
template<typename fp_t>
void tsrceKahanReal(vector<fp_t> &coefficients, vector<fp_t> &roots);
template<typename fp_t>
void tsrceKahan(vector<fp_t> &coefficients, vector<complex<fp_t>> &roots);

#endif //TSRCEKAHAN_TSRCEKAHAN_H
