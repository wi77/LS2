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

#include <math.h>
#include <float.h>

#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include "nlosp_em.h"

static float nlosp_mean   = 50.0F;
static float nlosp_sdev   = 15.0F;
static float nlosp_nlos_p = 0.1F;
static float nlosp_rate   = 2.0F;
static float nlosp_scale  = 100.0F;

#if defined(HAVE_POPT_H)
struct poptOption nlosp_arguments[] = {
  { "nlosp-nd-mean", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &nlosp_mean, 0,
    "mean value of the error", NULL },
  { "nlosp-nd-sdev", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &nlosp_sdev, 0,
    "standard deviation of the error", NULL },
  { "nlosp-prob", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &nlosp_nlos_p, 0,
    "Probability of an NLOS error", NULL },
  { "nlosp-rate", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &nlosp_rate, 0,
    "Rate of exponential NLOS error", NULL },
  { "nlosp-scale", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
    &nlosp_scale, 0,
    "Scale of exponential NLOS error", NULL },
  POPT_TABLEEND
};
#endif


void
nlosp_setup(const vector2 *anchors __attribute__((__unused__)),
            size_t nanchors __attribute__((__unused__)))
{
    // Do nothing.
}

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
        VECTOR rn = gaussrand(seed, nlosp_mean, nlosp_sdev);
        VECTOR help = VECTOR_LT(rnd(seed), VECTOR_BROADCASTF(nlosp_nlos_p));
        VECTOR nlos = exp_rand(seed, VECTOR_BROADCASTF(nlosp_rate)) *
                          VECTOR_BROADCASTF(nlosp_scale);
	rn += VECTOR_AND(nlos, help);  // Mask out members without nlos error
        result[k] = distances[k] + rn;
    }
}

