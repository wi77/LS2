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

#include "bahillo_em.h"

/* @error_model_name: NLOS (Bahillo) */


// Errormodels have to include all utils themselves
#include "../util/util_random.c"



extern void __attribute__((__nonnull__))
ls2_init_bahillo_arguments(ls2_bahillo_arguments *arguments)
{
        arguments->mean = 5.0;
        arguments->sdev = 2.5;
        arguments->ua = 2.38419e-07;
        arguments->ub = 3.0;
}



static ls2_bahillo_arguments bahillo_arguments;



static GOptionEntry bahillo_parameters[] = {
  { "bahillo-nd-mean", 0, 0, G_OPTION_ARG_DOUBLE,
    &bahillo_arguments.mean,
    "mean value of the error", NULL },
  { "bahillo-nd-sdev", 0, 0, G_OPTION_ARG_DOUBLE,
    &bahillo_arguments.sdev,
    "standard deviation of the error", NULL },
  { "bahillo-rate-low", 0, 0, G_OPTION_ARG_DOUBLE,
    &bahillo_arguments.ua,
    "Lower value of random NLOS rate", NULL },
  { "bahillo-rate-high", 0, 0, G_OPTION_ARG_DOUBLE,
    &bahillo_arguments.ub,
    "Upper value of random NLOS rate", NULL },
  { NULL }
};



void __attribute__((__nonnull__))
ls2_add_bahillo_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("bahillo",
                                "Parameters to the error model of Bahillo",
                                "Parameters to the error model of Bahillo",
                                NULL, NULL);
     g_option_group_add_entries(group, bahillo_parameters);
     g_option_context_add_group(context, group);
}



void
bahillo_setup(const vector2 *anchors __attribute__((__unused__)),
              size_t nanchors __attribute__((__unused__)))
{
    /* Do nothing. */
}



/* NLOS Error model proposed by Bahillo. 
 */
static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
bahillo_error(__m128i *restrict seed,
              const size_t no_anchors,
              const VECTOR *restrict distances,
              const VECTOR *restrict vx __attribute__((__unused__)),
              const VECTOR *restrict vy __attribute__((__unused__)),
              const VECTOR tagx __attribute__((__unused__)),
              const VECTOR tagy __attribute__((__unused__)),
              VECTOR *restrict result)
{
    static const VECTOR scale = VECTOR_CONST_BROADCAST(50.0f/6.8f);
    static const VECTOR threshold = VECTOR_CONST_BROADCAST(1000.0f);

    for (size_t k = 0; k < no_anchors ; k++) {
        VECTOR error;
	do {
            VECTOR noise, nlos, rate;
	    // Compute the measurement noise.
            noise = gaussrand(seed, (float) bahillo_arguments.mean, (float) bahillo_arguments.sdev);
	    // Compute the non-line-of-sight error.
            rate = rnd(seed) * VECTOR_BROADCASTF((float) (bahillo_arguments.ub - bahillo_arguments.ua)) +
                    VECTOR_BROADCASTF((float) bahillo_arguments.ua);
            nlos = exp_rand(seed, rate);
            error = (noise + nlos) * scale;
        } while(!VECTOR_TEST_ALL_ONES(VECTOR_LE(error, threshold)));

        // In Bahillo's paper, the nlos error uses a rate from (0, 3].
        // The expectation of the exponential is estimated to be 5.1.
        // Since we want an expected value of 50, we multiply the nlos
        // accordingly.
        result[k] = distances[k] + error;
    }
}

