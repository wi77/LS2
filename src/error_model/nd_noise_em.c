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

/* @error_model_name: Normal error */


#ifdef HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#include <glib.h>

#include "nd_noise_em.h"

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

static double nd_noise_mean = 50.0F;
static double nd_noise_sdev = 25.0F;

static GOptionEntry nd_noise_arguments[] = {
        { "nd-mean", 0, 0, G_OPTION_ARG_DOUBLE, &nd_noise_mean,
          "mean value of the error", NULL },
        { "nd-sdev", 0, 0, G_OPTION_ARG_DOUBLE, &nd_noise_sdev,
          "standard deviation of the error", NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_nd_noise_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("nd-noise",
                                "Parameters to the Gaussian noise error model",
                                "Parameters to the Gaussian noise error model",
                                NULL, NULL);
     g_option_group_add_entries(group, nd_noise_arguments);
     g_option_context_add_group(context, group);
}


void
nd_noise_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
  // No setup needed.
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
nd_noise_error(__m128i *restrict seed,
               const size_t anchors, 
               const VECTOR *restrict distances,
               const VECTOR *restrict vx __attribute__((__unused__)),
               const VECTOR *restrict vy __attribute__((__unused__)),
               const VECTOR tagx __attribute__((__unused__)),
               const VECTOR tagy __attribute__((__unused__)),
               VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
            VECTOR rn = gaussrand(seed, (float) nd_noise_mean,
                                  (float) nd_noise_sdev);
        //if negative Values are not acceptable, use this (about 10% slower, half a sec)
        /*
        while(rn[0]<0||rn[1]<0||rn[2]<0||rn[3]<0){
        
        rn = gaussrand(seed);
        
        }*/
      	result[k] = distances[k] + rn;
    }
}

