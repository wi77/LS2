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

/* @error_model_name: Gamma error (buggy) */

#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include "gamma_noise_em.h"

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

/* Mean = shape / rate
 * Variance = shape / (rate * rate)
 *
 * The defaults result in mean = 50.0f and sigma = 28.867
 */
static float gamma_shape = 3.0f;
static float gamma_rate = 3.0f / 50.0f;   // mean = shape / rate
static float gamma_offset = 0.0f;


#if defined(HAVE_POPT_H)
struct poptOption gamma_noise_arguments[] = {
    { "gamma-shape", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &gamma_shape, 0,
      "shape of the gamma distribution", NULL },
    { "gamma-rate", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &gamma_rate, 0,
      "rate of the gamma distribution", NULL },
    { "gamma-offset", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &gamma_offset, 0,
      "additive offset to the gamma distribution", NULL },
    POPT_TABLEEND
};
#endif


void
gamma_noise_setup(const vector2 *anchors __attribute__((__unused__)),
                  size_t nanchors __attribute__((__unused__)))
{
  // No setup needed.
}

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
gamma_noise_error(__m128i *restrict seed,
                  const size_t anchors,
                  const VECTOR *restrict distances,
                  const VECTOR *restrict vx __attribute__((__unused__)),
                  const VECTOR *restrict vy __attribute__((__unused__)),
                  const VECTOR tagx __attribute__((__unused__)),
                  const VECTOR tagy __attribute__((__unused__)),
                  VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
        VECTOR x = VECTOR_BROADCASTF(1.0F);
        float alpha;
	for (alpha = gamma_shape; alpha >= 1.0F; alpha -= 1.0F) {
            x *= rnd(seed);
        }
        if (alpha > 0.0F) {
            assert(alpha < 1.0F);
            // Ahrens-Dieter acceptance-rejection method
            const VECTOR v0 =
                VECTOR_BROADCASTF(((float) M_E) / ((float) M_E + alpha));
            VECTOR xi, eta, mask;
            do {
                // TODO: implement a way to keep accepted random values.
                // Here, we try until all numbers are accepted.
                VECTOR xi0, xi1, eta0, eta1;
                VECTOR V0 = rnd(seed), V1 = rnd(seed), V2 = rnd(seed);
                VECTOR m = VECTOR_LE(V0, v0);
                xi0 = VECTOR_POW(V1, one / alpha);
                eta0 = V2 * VECTOR_POW(xi0, VECTOR_BROADCASTF(alpha - 1.0F));
                xi1 = VECTOR_BROADCASTF(1.0F) - VECTOR_LOG(V1);
                eta1 = V2 * VECTOR_EXP(-xi1);
                xi = VECTOR_OR(VECTOR_AND(xi0, m), VECTOR_ANDNOT(xi1, m));
                eta = VECTOR_OR(VECTOR_AND(eta0, m), VECTOR_ANDNOT(eta1, m)); 
                mask = VECTOR_GT(eta, VECTOR_POW(xi, alpha - VECTOR_BROADCASTF(1.0F)) * VECTOR_EXP(-xi));
            } while (VECTOR_TEST_ALL_ONES(VECTOR_NE(mask, VECTOR_ZERO())));
            x *= xi;
        }
        x = (VECTOR_LOG(x) / VECTOR_BROADCASTF(-gamma_rate)) -
              VECTOR_BROADCASTF(gamma_offset);

      	result[k] = distances[k] + x;
    }
}
