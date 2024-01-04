//
// Created by elturpin on 10/12/2020.
//

#include "RandService.h"

#include <cassert>
#include <cmath>

__device__ int RandService::random_roulette(double* probability_array, int nb_elements, uint idx, Phase phase) {
    auto rand_number = gen_double(idx, phase);
    int selection = -1;
    do
    {
        assert(selection < nb_elements - 1);
        rand_number -= probability_array[++selection];
    } while (rand_number > 0);

    return selection;
}

static __constant__ double cof[6] = {76.18009172947146,
                                     -86.50532032941677,
                                     24.01409824083091,
                                     -1.231739572450155,
                                     0.1208650973866179e-2,
                                     -0.5395239384953e-5};

// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt).
// When the argument z is an integer, the gamma function is just the familiar factorial
// function, but offset by one, n! = gamma(n + 1).
__device__ static double gammln(double X)
{
    double x, y, tmp, ser;

    y = x = X;
    tmp = x + 5.5;
    tmp -= (x+0.5) * log(tmp);
    ser = 1.000000000190015;

    for (int8_t j = 0 ; j <= 5 ; j++)
    {
        ser += cof[j] / ++y;
    }

    return -tmp + log(2.5066282746310005 * ser / x);
}

__device__
int32_t RandService::binomial_random(int32_t nb_drawings, double prob, uint idx, Phase phase) {
    int32_t nb_success;

    // The binomial distribution is invariant under changing
    // ProbSuccess to 1-ProbSuccess, if we also change the answer to
    // NbTrials minus itself; we ll remember to do this below.
    double p;
    if (prob <= 0.5) p = prob;
    else p = 1.0 - prob;

    // mean of the deviate to be produced
    double mean = nb_drawings * p;


    if (nb_drawings < 25)
        // Use the direct method while NbTrials is not too large.
        // This can require up to 25 calls to the uniform random.
    {
        nb_success = 0;
        for (int32_t j = 1 ; j <= nb_drawings ; j++)
        {
            if (gen_double(idx, phase) < p) nb_success++;
        }
    }
    else if (mean < 1.0)
        // If fewer than one event is expected out of 25 or more trials,
        // then the distribution is quite accurately Poisson. Use direct Poisson method.
    {
        double g = exp(-mean);
        double t = 1.0;
        int32_t j;
        for (j = 0; j <= nb_drawings ; j++)
        {
            t = t * gen_double(idx, phase);
            if (t < g) break;
        }

        if (j <= nb_drawings) nb_success = j;
        else nb_success = nb_drawings;
    }

    else
        // Use the rejection method.
    {
        double en     = nb_drawings;
        double oldg   = gammln(en + 1.0);
        double pc     = 1.0 - p;
        double plog   = log(p);
        double pclog  = log(pc);

        // rejection method with a Lorentzian comparison function.
        double sq = sqrt(2.0 * mean * pc);
        double angle, y, em, t;
        do
        {
            do
            {
                angle = M_PI * gen_double(idx, phase);
                y = tan(angle);
                em = sq*y + mean;
            } while (em < 0.0 || em >= (en + 1.0)); // Reject.

            em = floor(em); // Trick for integer-valued distribution.
            t = 1.2 * sq * (1.0 + y*y)
                * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);

        } while (gen_double(idx, phase) > t); // Reject. This happens about 1.5 times per deviate, on average.

        nb_success = (int32_t) rint(em);
    }


    // Undo the symmetry transformation.
    if (p != prob) nb_success = nb_drawings - nb_success;

    return nb_success;
}
