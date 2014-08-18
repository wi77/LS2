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

#include <glib.h>
#include "erlang_noise_em.h"

static int erlang_shape = 3;
static double erlang_rate = 0.5228819579;
static double erlang_offset = 3.31060119642765;
static double erlang_scale = 50.0 / 2.85;

static GOptionEntry erlang_noise_arguments[] = {
    { "erlang-shape", 0, 0, G_OPTION_ARG_INT, &erlang_shape,
      "shape of the erlang distribution", NULL },
    { "erlang-rate", 0, 0, G_OPTION_ARG_DOUBLE, &erlang_rate,
      "rate of the erlang distribution", NULL },
    { "erlang-offset", 0, 0, G_OPTION_ARG_DOUBLE, &erlang_offset,
      "additive offset to the erlang distribution", NULL },
    { "erlang-scale", 0, 0, G_OPTION_ARG_DOUBLE, &erlang_scale,
      "multiplier to scale the erlang distribution", NULL },
    { NULL }
};


void __attribute__((__nonnull__))
ls2_add_erlang_noise_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("erlang-noise",
                                "Parameters to the Erlang noise error model",
                                "Parameters to the Erlang noise error model",
                                NULL, NULL);
     g_option_group_add_entries(group, erlang_noise_arguments);
     g_option_context_add_group(context, group);
}


void
erlang_noise_setup(const vector2 *anchors __attribute__((__unused__)),
                   size_t nanchors __attribute__((__unused__)))
{
  // No setup needed.
}


// Errormodels have to include all utils themselves
#if defined(STAND_ALONE)
#  include "../util/util_random.c"
#endif


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,3,8)))
erlang_noise_error(__m128i *restrict seed,
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
	for (int i = 0; i < erlang_shape; i++)
	    x *= rnd(seed);
        x = VECTOR_LOG(x) / VECTOR_BROADCASTF((float) -erlang_rate);
        x -= VECTOR_BROADCASTF((float) erlang_offset);

	// HACK: Scale it to an expected value of 50.
        x *= VECTOR_BROADCASTF((float) erlang_scale);
      	result[k] = distances[k] + x;
    }
}
