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

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

#include <assert.h>
#include <math.h>
#include <float.h>

#include "nlosp_em.h"

#ifndef ERROR
#  define ERROR 100
#endif

#ifndef NSUM
#  define NSUM 25
#endif

#ifndef Distribution
#  define Distribution 2
#endif

#ifndef NLOS
#  define NLOS 0.1F
#endif

VECTOR r_error;

void
nlosp_setup(const vector2 *anchors __attribute__((__unused__)),
            size_t nanchors __attribute__((__unused__)))
{
    r_error = VECTOR_BROADCASTF(ERROR);    
}

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
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
        VECTOR rn = gaussrand(seed, 50, 15);
        VECTOR help = VECTOR_LT(rnd(seed), VECTOR_BROADCASTF(NLOS));
        VECTOR nlos = exp_rand(seed, VECTOR_BROADCASTF(Distribution)) *
                      VECTOR_BROADCASTF(100);
	rn += VECTOR_AND(nlos, help);  // Mask out members without nlos error
        result[k] = distances[k] + rn;
    }
}

