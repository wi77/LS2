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

#include <glib.h>

#include "ab_nlos_em.h"

/* @error_model_name: Anchor Based NLOS error */
// Errormodels have to include all utils themselves
#include "../util/util_random.c"

static double ab_nlos_mean   = 50.0F;
static double ab_nlos_sdev   = 15.0F;
static int ab_nlos_count    = 3;
static double ab_nlos_rate   = 2.0F;
static double ab_nlos_scale  = 100.0F;
static VECTOR ab_nlos_oscale;
static int ab_nlos_norm = 1;

static GOptionEntry ab_nlos_arguments[] = {
  { "ab_nlos-nd-mean", 0, 0, G_OPTION_ARG_DOUBLE,
    &ab_nlos_mean,
    "mean value of the error", NULL },
  { "ab_nlos-nd-sdev", 0, 0, G_OPTION_ARG_DOUBLE,
    &ab_nlos_sdev,
    "standard deviation of the error", NULL },
  { "ab_nlos-count", 0, 0, G_OPTION_ARG_INT,
    &ab_nlos_count,
    "First n anchors with NLOS ERROR", NULL },
  { "ab_nlos-norm", 0, 0, G_OPTION_ARG_INT,
    &ab_nlos_norm,
    "Normalize error?", NULL },
  { "ab_nlos-rate", 0, 0, G_OPTION_ARG_DOUBLE,
    &ab_nlos_rate,
    "Rate of exponential NLOS error", NULL },
  { "ab_nlos-scale", 0, 0, G_OPTION_ARG_DOUBLE,
    &ab_nlos_scale,
    "Scale of exponential NLOS error", NULL },
  { NULL }
};


void __attribute__((__nonnull__))
ls2_add_ab_nlos_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("ab-nlos",
                                "Parameters to the AB NLOS error model",
                                "Parameters to the AB NLOS error model",
                                NULL, NULL);
     g_option_group_add_entries(group, ab_nlos_arguments);
     g_option_context_add_group(context, group);
}


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
                rn = gaussrand(seed, (float) ab_nlos_mean, (float) ab_nlos_sdev);
                VECTOR nlos = exp_rand(seed, VECTOR_BROADCASTF((float) ab_nlos_rate)) * VECTOR_BROADCASTF((float) ab_nlos_scale);
	        rn += nlos;               
        } else {
                rn = gaussrand(seed, (float) ab_nlos_mean, (float) ab_nlos_sdev);
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
    if (!ab_nlos_norm) return;
    
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
