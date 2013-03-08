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

/* @error_model_name: Constant error */

#include "const_em.h"

static float const_error_value = 50.0F;

#ifdef HAVE_POPT_H
struct poptOption const_arguments[] = {
        { "const-error", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &const_error_value, 0,
          "constant error", NULL },
        POPT_TABLEEND
};
#endif


static VECTOR error_v;


void
const_setup(const vector2 *anchors __attribute__((__unused__)),
            size_t nanchors __attribute__((__unused__)))
{
    error_v = VECTOR_BROADCASTF(const_error_value);
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
const_error(__m128i *restrict seed __attribute__((__unused__)),
            const size_t anchors,
            const VECTOR *restrict distances,
            const VECTOR *restrict vx __attribute__((__unused__)),
            const VECTOR *restrict vy __attribute__((__unused__)),
            const VECTOR tagx __attribute__((__unused__)),
            const VECTOR tagy __attribute__((__unused__)),
            VECTOR *restrict result)
{
    for (size_t k = 0; k < anchors ; k++) {
      	result[k] = distances[k] + error_v;
    }
}
