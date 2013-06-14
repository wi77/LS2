/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
  Stefan Adler, Malte Rohde, Jonathan Gunthermann
 
  LS² is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  LS² is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with LS².  If not, see <http://www.gnu.org/licenses/>.

 */

/* @error_model_name: Anchor Based NLOS error */
// Errormodels have to include all utils themselves
#include "../util/util_random.c"
#include "ab_nlos_em.h"

static float ab_nlos_mean   = 50.0F;
static float ab_nlos_sdev   = 15.0F;
static int ab_nlos_count    = 3;
static float ab_nlos_rate   = 2.0F;
static float ab_nlos_scale  = 100.0F;
static VECTOR ab_nlos_oscale;

#if defined(HAVE_POPT_H)
struct poptOption ab_nlos_arguments[] = {
  { "ab_nlos-nd-mean", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &ab_nlos_mean, 0,
    "mean value of the error", NULL },
  { "ab_nlos-nd-sdev", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &ab_nlos_sdev, 0,
    "standard deviation of the error", NULL },
  { "ab_nlos-count", 0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
    &ab_nlos_count, 0,
    "First n anchors with NLOS ERROR", NULL },
  { "ab_nlos-rate", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &ab_nlos_rate, 0,
    "Rate of exponential NLOS error", NULL },
  { "ab_nlos-scale", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &ab_nlos_scale, 0,
    "Scale of exponential NLOS error", NULL },
  POPT_TABLEEND
};
#endif

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
ab_nlos_error(__m128i *restrict seed,
            const size_t anchors,
            const VECTOR *restrict distances,
            const VECTOR *restrict vx __attribute__((__unused__)),
            const VECTOR *restrict vy __attribute__((__unused__)),
            const VECTOR tagx __attribute__((__unused__)),
            const VECTOR tagy __attribute__((__unused__)),
            VECTOR *restrict result)
{
    VECTOR rn;
    for (size_t k=0; k < anchors ; k++) {
        if (k < (size_t)ab_nlos_count) {
            rn = gaussrand(seed, ab_nlos_mean, ab_nlos_sdev);
            VECTOR nlos = exp_rand(seed, VECTOR_BROADCASTF(ab_nlos_rate)) * VECTOR_BROADCASTF(ab_nlos_scale);
	        rn += nlos;               
        } else {
            rn = gaussrand(seed, ab_nlos_mean, ab_nlos_sdev);
	    }
        result[k] = distances[k] + (rn*ab_nlos_oscale);
    }
}

void
ab_nlos_setup(const vector2 *anchors __attribute__((__unused__)), size_t nanchors)
{
    #define TESTRUNS 10000.0f
    const VECTOR d;
    VECTOR test[nanchors];
    ab_nlos_oscale = one;
    for (int i = 0; i < (int)nanchors; i++)
        test[i]=zero;
    unsigned int _seed = (unsigned int)time(NULL);
    double mean=0.0;
    __m128i seed;

    int seed0 = rand_r(&_seed);
    int seed1 = rand_r(&_seed);
    int seed2 = rand_r(&_seed);
    int seed3 = rand_r(&_seed);
    seed = _mm_set_epi32(seed0, seed1, seed2, seed3);

    for (int i=0; i < TESTRUNS; i++){
        ab_nlos_error(&seed,nanchors,test,&d,&d,d,d,test);
        for (int a=0; a<(int)nanchors; a++) {
            test[a] /= TESTRUNS;
            for (int ii=0; ii < VECTOR_OPS; ii++) {
                mean += (test[a][ii] / (double)VECTOR_OPS) / (double)nanchors;
            }
            test[a]=zero;
        }
    }
    float scale = (float)(50.0 / mean);
    ab_nlos_oscale = VECTOR_BROADCASTF(scale);    
}
