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

#include "bahillo_em.h"

/* @error_model_name: NLOS error (Normal + Exponential) */

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

void
bahillo_setup(const vector2 *anchors __attribute__((__unused__)),
              size_t nanchors __attribute__((__unused__)))
{
    /* Do nothing. */
}


static const float mean = 5.0f;
static const float sdev = 2.5f;

/* The next two parameters specify the range of rates. Errors can be
 * huge if the rate is close to zero, in order to avoid those huge
 * values.
 */
static float bahillo_ua = 2.38419e-07f;
static float bahillo_ub = 3.0f;


/* NLOS Error model proposed by Bahillo. 
 */
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
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
            noise = gaussrand(seed, mean, sdev);
	    // Compute the non-line-of-sight error.
            rate = rnd(seed) * VECTOR_BROADCASTF(bahillo_ub - bahillo_ua) +
                   VECTOR_BROADCASTF(bahillo_ua);
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

