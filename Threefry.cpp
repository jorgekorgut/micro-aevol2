#include "Threefry.h"

#include <cmath>

void Threefry::save(gzFile backup_file) const {
  Threefry123::key_type::value_type seed = seed_[1];
  gzwrite(backup_file, &seed, sizeof(seed));
  gzwrite(backup_file, counters_.data(), counters_.size() * sizeof(counters_[0]));
}

Threefry::Threefry(int X, int Y, gzFile backup_file)
        : counters_(X * Y * NPHASES, 0), X_(X), Y_(Y), N_(X * Y) {

  Threefry123::key_type::value_type seed;
  gzread(backup_file, &seed, sizeof(seed));
  seed_[0] = 0;
  seed_[1] = seed;

  crt_value_type tmp_counters[counters_.size()];
  gzread(backup_file, tmp_counters, counters_.size() * sizeof(tmp_counters[0]));

  counters_ = std::vector<crt_value_type>(tmp_counters, tmp_counters + counters_.size());
}

int32_t Threefry::Gen::roulette_random(double* probs, int32_t nb_elts, bool verbose )
{
    //cloned_probs.resize(nb_elts);
    //for (int i = 0; i < nb_elts; i++) cloned_probs[i] = probs[i];

  double pick_one = 0.0;

  while (pick_one == 0.0)
  {
    pick_one = random();
    //pickones.push_back(pick_one);
    //if (verbose) printf("pick one : %f\n",pick_one);
  }

  int32_t found_org = 0;

  pick_one -= probs[0];
  while (pick_one > 0)
  {
    assert(found_org<nb_elts-1);
    //pickones3.push_back(probs[found_org+1]);

    pick_one -= probs[++found_org];
    //pickones2.push_back(pick_one);
  }
  return found_org;
}

// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt).
// When the argument z is an integer, the gamma function is just the familiar factorial
// function, but offset by one, n! = gamma(n + 1).
static double gammln(double X)
{
  double x, y, tmp, ser;
  static double cof[6] = {  76.18009172947146,
                            -86.50532032941677,
                            24.01409824083091,
                            -1.231739572450155,
                            0.1208650973866179e-2,
                            -0.5395239384953e-5 };

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

/*!
  Binomial drawing of parameter (nb_drawings, prob).

  Number of successes out of nb_drawings trials each of probability prob.
 */
int32_t Threefry::Gen::binomial_random(int32_t nb_drawings, double prob)
{
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
      if (random() < p) nb_success++;
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
      t = t * random();
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
        angle = M_PI * random();
        y = tan(angle);
        em = sq*y + mean;
      } while (em < 0.0 || em >= (en + 1.0)); // Reject.

      em = floor(em); // Trick for integer-valued distribution.
      t = 1.2 * sq * (1.0 + y*y)
              * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);

    } while (random() > t); // Reject. This happens about 1.5 times per deviate, on average.

    nb_success = (int32_t) rint(em);
  }


  // Undo the symmetry transformation.
  if (p != prob) nb_success = nb_drawings - nb_success;

  return nb_success;
}
