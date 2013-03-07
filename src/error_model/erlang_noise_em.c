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

/* @error_model_name: Erlang error */

#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include "erlang_noise_em.h"

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

static int erlang_shape = 3;
static float erlang_rate = 0.5228819579F;
static float erlang_offset = 3.31060119642765F;
static float erlang_scale = 50.0f / 2.85f;

#if defined(HAVE_POPT_H)
struct poptOption erlang_noise_arguments[] = {
    { "erlang-shape", 0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
      &erlang_shape, 0,
      "shape of the erlang distribution", NULL },
    { "erlang-rate", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &erlang_rate, 0,
      "rate of the erlang distribution", NULL },
    { "erlang-offset", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &erlang_offset, 0,
      "additive offset to the erlang distribution", NULL },
    { "erlang-scale", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
      &erlang_scale, 0,
      "multiplier to scale the erlang distribution", NULL },
    POPT_TABLEEND
};
#endif

void
erlang_noise_setup(const vector2 *anchors __attribute__((__unused__)),
                   size_t nanchors __attribute__((__unused__)))
{
  // No setup needed.
}

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
erlang_noise_error(__m128i* seed,
                   const size_t anchors,
                   const VECTOR *restrict distances __attribute__((__unused__)),
                   const VECTOR *restrict vx __attribute__((__unused__)),
                   const VECTOR *restrict vy __attribute__((__unused__)),
                   const VECTOR tagx __attribute__((__unused__)),
                   const VECTOR tagy __attribute__((__unused__)),
                   VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
	VECTOR x = VECTOR_BROADCASTF(1.0F);
	for (int i = 0; i < erlang_shape; i++)
	    x *= rnd(seed);
        x = VECTOR_LOG(x) / VECTOR_BROADCASTF(-erlang_rate);
        x -= VECTOR_BROADCASTF(erlang_offset);

	// HACK: Scale it to an expected value of 50.
        x *= VECTOR_BROADCASTF(erlang_scale);
      	result[k] = distances[k] + x;
    }
}
