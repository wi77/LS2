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

/* @error_model_name: Weibull noise */


#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include <float.h>

#include "weibull_em.h"

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

/* Mean value is scale * sqrtf(M_PI / 2.0f). */
static float weibull_scale = 55.571f;
static float weibull_shape = 3.5f;

static VECTOR weibull_scale_vector;
static VECTOR weibull_shape_vector;

#if defined(HAVE_POPT_H)
struct poptOption weibull_arguments[] = {
        { "weibull-scale", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &weibull_scale, 0,
          "scale parameter of the Weibull noise", NULL },
        { "weibull-shape", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &weibull_shape, 0,
          "shape parameter of the Weibull noise", NULL },
        POPT_TABLEEND
};
#endif

void
weibull_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
    VECTOR t = VECTOR_CONST_BROADCAST(weibull_scale);
    weibull_scale_vector = t;
    VECTOR u = VECTOR_CONST_BROADCAST(weibull_shape);
    weibull_shape_vector = u;
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
weibull_error(__m128i *restrict seed,
               const size_t anchors, 
               const VECTOR *restrict distances,
               const VECTOR *restrict vx __attribute__((__unused__)),
               const VECTOR *restrict vy __attribute__((__unused__)),
               const VECTOR tagx __attribute__((__unused__)),
               const VECTOR tagy __attribute__((__unused__)),
               VECTOR *restrict result)
{
    static const VECTOR low = VECTOR_CONST_BROADCAST(FLT_EPSILON);
    static const VECTOR high = VECTOR_CONST_BROADCAST(1.0f - FLT_EPSILON);
    for (size_t k = 0; k < anchors ; k++) {
        VECTOR t = rnd(seed); // Uniform distributed in [0,1]
        t = VECTOR_MAX(t, low);
        t = VECTOR_MIN(t, high);
        t = one - t;
        t = VECTOR_LOG(t);
        t = zero - t;
        t = VECTOR_POW(t, one / weibull_shape_vector);
        t *= weibull_scale_vector;
      	result[k] = distances[k] + t;
    }
}

