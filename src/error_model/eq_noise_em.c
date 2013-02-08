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

#include "eq_noise_em.h"

float eq_error_min = 0.0F;
float eq_error_max = 100.0F;

#ifdef HAVE_POPT_H
struct poptOption eq_arguments[] = {
        { "eq-error-min", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &eq_error_min, 0,
          "minimum value of error", NULL },
        { "eq-error-max", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &eq_error_max, 0,
          "maximum value of error", NULL },
        POPT_TABLEEND
};
#endif

// Errormodels have to include all utils themselves
#if defined(STAND_ALONE)
#  include "../util/util_random.c"
#endif

VECTOR eq_error_min_v;
VECTOR eq_error_rng_v;

void
eq_noise_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
    eq_error_min_v = VECTOR_BROADCASTF(eq_error_min);
    eq_error_rng_v = VECTOR_BROADCASTF(eq_error_max - eq_error_min);
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
eq_noise_error(__m128i* seed,
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
