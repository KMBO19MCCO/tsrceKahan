#include "excerpt.h"

template<typename fp_t>
int generate_polynomial(
        unsigned P, // polynomial degree
        unsigned N_pairs_of_complex_roots, // how many pairs of complex conjugate roots to introduce
        unsigned N_clustered_roots, // how many clustered roots to introduce; all the clustered roots are real
        unsigned N_multiple_roots, // how many multiple roots to introduce; all multiple roots are real
        fp_t max_distance_between_clustered_roots, // maximal distance between the closest of the clustered roots
        fp_t root_sweep_low,
        fp_t root_sweep_high, // low and high boundaries of real roots; imaginary parts of complex conjugate roots are in the same range
        std::vector<fp_t> &roots, // storage where to put the roots; size should exceed P-1
        std::vector<fp_t> &coefficients) // storage where to put the coefficients; size should exceed P
{
    int n_simple_roots = P - 2 * N_pairs_of_complex_roots - N_clustered_roots - N_multiple_roots;
    assert(N_clustered_roots != 1);
    assert(N_multiple_roots != 1);
    assert(n_simple_roots >= 0 && P <= 4);
    assert(root_sweep_high - root_sweep_low > 2 * P * max_distance_between_clustered_roots);

    std::random_device r;

    unsigned long long seed =
            std::chrono::system_clock::now().time_since_epoch().count() + r(); // counts milliseconds
    std::mt19937_64 rng(seed); // randomize seed from the clock
    std::uniform_real_distribution<fp_t> rnr(root_sweep_low,
                                             root_sweep_high); // uniform random data generator for single roots
    std::uniform_real_distribution<fp_t> rnc(static_cast<fp_t>(0.0L),
                                             max_distance_between_clustered_roots); // uniform random data generator for root clusters
    fp_t re, im, u, v = rnr(rng), root_mid_sweep = root_sweep_low + 0.5 * (root_sweep_high - root_sweep_low);

    coefficients[P] = static_cast<fp_t>(1.0L); // invariant
    switch (P) {
        case 0:
            coefficients[0] = rnr(rng);
            return 0;
        case 1:
            coefficients[0] = -(roots[0] = rnr(rng));
            return 1;
        case 2: {
            if (N_pairs_of_complex_roots == 1) // no real roots
            {
                re = rnr(rng);
                while ((im = rnr(rng)) == static_cast<fp_t>(0.0L)) {}
                coefficients[1] = re * static_cast<fp_t>(-2.0L);
//                coefficients[0] = pr_product_difference(re, re, -im, im); // re*re+im*im
                return 0;
            } else if (N_clustered_roots == 2) // 2 close but distinct roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re >= root_mid_sweep ? re - im : re + im);
            } else if (N_multiple_roots == 2) // double root counted as a single root
            {
                roots[1] = roots[0] = im = re = rnr(rng);
            } else // 2 distinct single roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnr(rng)) == re) {}
                roots[1] = im = rnr(rng);
            }
            coefficients[1] = -re - im;
            coefficients[0] = re * im;
            return 2; // return ((re!=im) ? 2 : 1);
        } // P=2
        case 3: {
            if (N_pairs_of_complex_roots == 1) // one real root
            {
                re = rnr(rng);
                while ((im = rnr(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[0] = u = rnr(rng);
//                im = pr_product_difference(re, re, -im, im); // re*re+im*im
                re *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by (x-u)
                coefficients[0] = -im * u;
                coefficients[2] = re - u;
                coefficients[1] = fma(-re, u, im); // im-re*u;
                return 1;
            } else if (N_clustered_roots == 3) // 3 clustered distinct roots
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                while ((u = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re > root_mid_sweep ? roots[2] = u = (re - im - u), re - im : roots[2] = u = (re + im +
                                                                                                               u), re +
                                                                                                                   im);
            } else if (N_clustered_roots == 2) // 2 clustered roots, 1 single root; all distinct
            {
                roots[0] = re = rnr(rng);
                while ((im = rnc(rng)) == static_cast<fp_t>(0.0L)) {}
                roots[1] = im = (re > root_mid_sweep ? re - im : re + im);
                do { roots[2] = u = rnr(rng); } while (u == re || u == roots[1]);
            } else if (N_multiple_roots == 3) // triple root counted as a single root
            {
                roots[2] = roots[1] = roots[0] = u = im = re = rnr(rng);
            } else if (N_multiple_roots == 2) // double root and 1 single root; totally 2 roots
            {
                roots[1] = roots[0] = im = re = rnr(rng);
                while ((roots[2] = u = rnr(rng)) == re) {}
            } else // 3 distinct single roots
            {
                roots[0] = re = rnr(rng);
                while ((roots[1] = im = rnr(rng)) == re) {}
                do { roots[2] = u = rnr(rng); } while (u == re || u == im);
            }
            coefficients[2] = -re - im - u;
            coefficients[0] = -re * im * u;
//            v = pr_product_difference(re, im, -re, u);
            coefficients[1] = fma(im, u, v); // re*im+re*u+im*u=im*u+(re*im-(-re*u));
            // if (re!=im && im!=u && u!=re) return 3;
            // if (re==im && im==u) return 1;
            // return 2;
            return 3;
        } // P=3
        case 4: // DEN DEBUG: check it carefully
        {
            if (N_pairs_of_complex_roots == 2) // no real roots
            {
                re = rnr(rng);
                im = rnr(rng);
                im = re * re + im * im;
                re *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + re*x + im)
                u = rnr(rng);
                v = rnr(rng);
                v = u * u + v * v;
                u *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + u*x + v)
                // multiply both irreducible quadratics
                coefficients[0] = im * v;
                coefficients[1] = re * v + im * u;
                coefficients[2] = im + re * u + v;
                coefficients[3] = re + u;
                return 0;
            } else if (N_pairs_of_complex_roots == 1) // two real roots
            {
                re = rnr(rng);
                im = rnr(rng);
                roots[0] = u = rnr(rng);
                im = re * re + im * im;
                re *= static_cast<fp_t>(-2.0L); // irreducible quadratic polynomial is (x^2 + re*x + im); multiply it by (x-u)

                if (N_clustered_roots == 2) // 2 clustered roots
                {
                    roots[0] = u = rnr(rng);
                    v = rnc(rng);
                    roots[1] = v = (u > root_mid_sweep ? u - v : u + v);
                } else if (N_multiple_roots == 2) // 2 multiple roots
                { roots[1] = roots[0] = u = v = rnr(rng); }
                else // 2 distinct roots
                {
                    roots[0] = u = rnr(rng);
                    roots[1] = v = rnr(rng);
                }
                root_mid_sweep = -u - v;
                v *= u;
                u = root_mid_sweep; // two-root quadratic polynomial is (x^2 + u*x + v)
                // multiply irreducible and reducible quadratics
                coefficients[0] = im * v;
                coefficients[1] = re * v + im * u;
                coefficients[2] = im + re * u + v;
                coefficients[3] = re + u;
                return 2;
            } else if (N_clustered_roots == 4) // 4 clustered roots
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                u = rnc(rng);
                v = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? roots[3] = v = (re - im - u - v), roots[2] = u = (re - im - u),
                        re - im :
                        roots[3] = v = (re + im + u + v), roots[2] = u = (re + im + u), re + im);
            } else if (N_clustered_roots == 3) // 3 clustered roots and 1 single root
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                u = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? roots[2] = u = (re - im - u), re - im : roots[2] = u = (re + im +
                                                                                                               u), re +
                                                                                                                   im);
                roots[3] = v = rnc(rng); // 1 single root
            } else if (N_clustered_roots == 2) // 2 clustered roots
            {
                roots[0] = re = rnr(rng);
                im = rnc(rng);
                roots[1] = im = (re > root_mid_sweep ? re - im : re + im);
                if (N_multiple_roots == 2) // 2 multiple roots
                { roots[3] = roots[2] = v = u = rnr(rng); }
                else // 2 single roots
                {
                    roots[2] = u = rnr(rng);
                    roots[3] = v = rnr(rng);
                }
            } else if (N_multiple_roots == 4) // 4 multiple roots
            { roots[3] = roots[2] = roots[1] = roots[0] = v = u = im = re = rnr(rng); }
            else if (N_multiple_roots == 3) // 3 multiple roots and 1 single root
            {
                roots[2] = roots[1] = roots[0] = u = im = re = rnr(rng);
                roots[3] = v = rnr(rng);
            } else if (N_multiple_roots == 2) // 2 multiple roots and 2 single roots
            {
                roots[1] = roots[0] = im = re = rnr(rng);
                roots[2] = u = rnr(rng);
                roots[3] = v = rnr(rng);
            } else // 4 distinct single roots
            {
                roots[0] = re = rnr(rng);
                roots[1] = im = rnr(rng);
                roots[2] = u = rnr(rng);
                roots[3] = v = rnr(rng);
            }
            // compute coefficients from 4 roots: re, im, u, v
            root_mid_sweep = -re - im;
            im *= re;
            re = root_mid_sweep; // now we have the 1.st quadratic polynomial: x^2 + x*re + im
            root_mid_sweep = -u - v;
            v *= u;
            u = root_mid_sweep; // now we have the 2.nd quadratic polynomial: x^2 + x*u + v
            coefficients[0] = im * v;
            coefficients[1] = re * v + im * u;
            coefficients[2] = re * u + im + v;
            coefficients[3] = re + u;
            return 4;
        } // P=4
        default:
            return -1;
    } // switch (P)
    return -1; // unreachable, means a flaw in control here
}

template<typename fp_t>
int compare_roots(
        unsigned N_roots_to_check, // number of roots in roots_to_check
        unsigned N_roots_ground_truth,  // number of roots in roots_ground_truth
        std::vector<fp_t> &roots_to_check, // one should take into account only first (N_roots_to_check) rots
        std::vector<fp_t> &roots_ground_truth, // one should take into account only first (N_roots_ground_truth) rots
        fp_t &max_deviation) // here will be placed the greatest among the smallest deviations of the roots in
// (roots_to_check) and (roots_ground_truth)
{
    int rv = PR_NUMBERS_OF_ROOTS_EQUAL;
    fp_t deviation, deviation_min_for_this_root, deviation_max = static_cast<fp_t>(0.0L);
    if (N_roots_to_check < N_roots_ground_truth) rv = PR_AT_LEAST_ONE_ROOT_LOST;
    else if (N_roots_to_check > N_roots_ground_truth) rv = PR_AT_LEAST_ONE_ROOT_IS_FAKE;
    // find the largest distance between the closest pairs of roots: one - from ground truth, one - from found ones
    for (int i = 0; i < N_roots_ground_truth; ++i) {
        // find the closest found root to the given ground truth root
        deviation_min_for_this_root = std::numeric_limits<fp_t>::infinity();
        for (int j = 0; j < N_roots_to_check; ++j) {
            deviation = std::abs(roots_ground_truth[i] - roots_to_check[j]);
            deviation_min_for_this_root =
                    deviation < deviation_min_for_this_root ? deviation : deviation_min_for_this_root;
        }
        deviation_max = deviation_min_for_this_root > deviation_max ? deviation_min_for_this_root : deviation_max;
    }
    max_deviation = deviation_max;
    return rv;
}

template int generate_polynomial<float>(unsigned P, unsigned N_pairs_of_complex_roots, unsigned N_clustered_roots,
                                        unsigned N_multiple_roots, float max_distance_between_clustered_roots,
                                        float root_sweep_low, float root_sweep_high, std::vector<float> &roots,
                                        std::vector<float> &coefficients);

