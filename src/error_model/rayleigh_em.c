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

/* @error_model_name: Rayleigh noise */


#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#include <glib.h>

#include <float.h>

#include "rayleigh_em.h"

/* Mean value is scale * sqrtf(M_PI / 2.0f). */
extern void __attribute__((__nonnull__))
ls2_init_rayleigh_arguments(ls2_rayleigh_arguments *arguments)
{
        arguments->scale = 50.0 / sqrt(M_PI / 2.0);
}



static ls2_rayleigh_arguments rayleigh_arguments;

static GOptionEntry rayleigh_parameters[] = {
        { "rayleigh-scale", 0, 0, G_OPTION_ARG_DOUBLE,
          &rayleigh_arguments.scale, "scale parameter of the Rayleigh noise",
          NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_rayleigh_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("rayleigh",
                                "Parameters to the Rayleigh error model",
                                "Parameters to the Rayleigh error model",
                                NULL, NULL);
     g_option_group_add_entries(group, rayleigh_parameters);
     g_option_context_add_group(context, group);
}


// Errormodels have to include all utils themselves
#include "../util/util_random.c"

static VECTOR rayleigh_scale_vector;

void
rayleigh_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
        VECTOR t = VECTOR_CONST_BROADCAST((float) rayleigh_arguments.scale);
        rayleigh_scale_vector = t;
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
rayleigh_error(__m128i *restrict seed,
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
    static const VECTOR minus_two = VECTOR_CONST_BROADCAST(-2.0f);
    for (size_t k = 0; k < anchors ; k++) {
        VECTOR t = rnd(seed); // Uniform distributed in [0,1]
        t = VECTOR_MAX(t, low);
        t = VECTOR_MIN(t, high);
        t = VECTOR_LOG(t);
        t *= minus_two;
        t = VECTOR_SQRT(t);
        t *= rayleigh_scale_vector;
      	result[k] = distances[k] + t;
    }
}
