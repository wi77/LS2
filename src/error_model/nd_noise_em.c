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

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include "nd_noise_em.h"

// Errormodels have to include all utils themselves
#include "../util/util_random.c"

static float nd_noise_mean = 50.0F;
static float nd_noise_sdev = 25.0F;

#if defined(HAVE_POPT_H)
struct poptOption nd_noise_arguments[] = {
        { "nd-mean", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &nd_noise_mean, 0,
          "mean value of the error", NULL },
        { "nd-sdev", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &nd_noise_sdev, 0,
          "standard deviation of the error", NULL },
        POPT_TABLEEND
};
#endif


void
nd_noise_setup(const vector2 *anchors __attribute__((__unused__)),
               size_t nanchors __attribute__((__unused__)))
{
  // No setup needed.
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
nd_noise_error(__m128i* seed,
               const size_t anchors, 
               const VECTOR *restrict distances,
               const VECTOR *restrict vx __attribute__((__unused__)),
               const VECTOR *restrict vy __attribute__((__unused__)),
               const VECTOR tagx __attribute__((__unused__)),
               const VECTOR tagy __attribute__((__unused__)),
               VECTOR *restrict result)
{
    for (size_t k=0; k < anchors ; k++) {
        VECTOR rn = gaussrand(seed, nd_noise_mean, nd_noise_sdev);
        //if negative Values are not acceptable, use this (about 10% slower, half a sec)
        /*
        while(rn[0]<0||rn[1]<0||rn[2]<0||rn[3]<0){
        
        rn = gaussrand(seed);
        
        }*/
      	result[k] = distances[k] + rn;
    }
}

