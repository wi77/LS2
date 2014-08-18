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

/* @error_model_name: NLOS error (Normal + Exponential) */

#include <math.h>
#include <float.h>

#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#include <glib.h>

#include "nlosp_em.h"

static double nlosp_mean   = 50.0;
static double nlosp_sdev   = 15.0;
static double nlosp_nlos_p = 0.1;
static double nlosp_rate   = 2.0;
static double nlosp_scale  = 100.0;

static GOptionEntry nlosp_arguments[] = {
  { "nlosp-nd-mean", 0, 0, G_OPTION_ARG_DOUBLE, &nlosp_mean,
    "mean value of the error", NULL },
  { "nlosp-nd-sdev", 0, 0, G_OPTION_ARG_DOUBLE, &nlosp_sdev,
    "standard deviation of the error", NULL },
  { "nlosp-prob", 0, 0, G_OPTION_ARG_DOUBLE, &nlosp_nlos_p,
    "Probability of an NLOS error", NULL },
  { "nlosp-rate", 0, 0, G_OPTION_ARG_DOUBLE, &nlosp_rate,
    "Rate of exponential NLOS error", NULL },
  { "nlosp-scale", 0, 0, G_OPTION_ARG_DOUBLE, &nlosp_scale,
    "Scale of exponential NLOS error", NULL },
  { NULL }
};


void __attribute__((__nonnull__))
ls2_add_nlosp_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("nlosp",
                                "Parameters to the NLOSP error model",
                                "Parameters to the NLOSP error model",
                                NULL, NULL);
     g_option_group_add_entries(group, nlosp_arguments);
     g_option_context_add_group(context, group);
}


void
nlosp_setup(const vector2 *anchors __attribute__((__unused__)),
            size_t nanchors __attribute__((__unused__)))
{
    // Do nothing.
}


// Errormodels have to include all utils themselves
#if defined(STAND_ALONE)
#  include "../util/util_random.c"
#endif

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
nlosp_error(__m128i *restrict seed,
            const size_t anchors,
            const VECTOR *restrict distances,
            const VECTOR *restrict vx __attribute__((__unused__)),
            const VECTOR *restrict vy __attribute__((__unused__)),
            const VECTOR tagx __attribute__((__unused__)),
            const VECTOR tagy __attribute__((__unused__)),
            VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
            VECTOR rn = gaussrand(seed, (float) nlosp_mean, (float) nlosp_sdev);
        VECTOR help = VECTOR_LT(rnd(seed), VECTOR_BROADCASTF((float) nlosp_nlos_p));
        VECTOR nlos = exp_rand(seed, VECTOR_BROADCASTF((float) nlosp_rate)) *
                          VECTOR_BROADCASTF((float) nlosp_scale);
	rn += VECTOR_AND(nlos, help);  // Mask out members without nlos error
        result[k] = distances[k] + rn;
    }
}
