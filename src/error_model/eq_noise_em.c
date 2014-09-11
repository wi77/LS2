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

/* @error_model_name: Uniform error */

#include <glib.h>
#include "eq_noise_em.h"

double eq_error_min = 0.0F;
double eq_error_max = 100.0F;

static GOptionEntry eq_noise_arguments[] = {
        { "eq-error-min", 0, 0, G_OPTION_ARG_DOUBLE, &eq_error_min,
          "minimum value of error", NULL },
        { "eq-error-max", 0, 0, G_OPTION_ARG_DOUBLE, &eq_error_max,
          "maximum value of error", NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_eq_noise_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("eq-noise",
                                "Parameters to the uniform distributed noise error model",
                                "Parameters to the uniform distributed noise error model",
                                NULL, NULL);
     g_option_group_add_entries(group, eq_noise_arguments);
     g_option_context_add_group(context, group);
}


// Errormodels have to include all utils themselves
#include "../util/util_random.c"

VECTOR eq_error_min_v;
VECTOR eq_error_rng_v;

void
eq_noise_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
        eq_error_min_v = VECTOR_BROADCASTF((float) eq_error_min);
        eq_error_rng_v = VECTOR_BROADCASTF((float) (eq_error_max - eq_error_min));
}


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
eq_noise_error(__m128i *restrict seed,
               const size_t anchors,
	       const VECTOR *restrict  const distances,
	       const VECTOR __attribute__((unused)) vx[MAX_ANCHORS],
	       const VECTOR __attribute__((unused)) vy[MAX_ANCHORS],
	       const VECTOR __attribute__((unused)) tagx,
	       const VECTOR __attribute__((unused)) tagy,
	       VECTOR *restrict result)
{
    for (size_t k = 0; k < anchors ; k++) {
        VECTOR rn = rnd(seed);
	rn *= eq_error_rng_v;
	rn += eq_error_min_v;
        result[k] = distances[k] + rn;
    }
}
