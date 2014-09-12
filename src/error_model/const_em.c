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

#include <glib.h>

#include "const_em.h"




extern void __attribute__((__nonnull__))
ls2_init_const_arguments(ls2_const_arguments *arguments)
{
        arguments->value = 50.0;
}



static ls2_const_arguments const_arguments;

static GOptionEntry const_parameters[] = {
        { "const-error", 0, 0, G_OPTION_ARG_DOUBLE, &const_arguments.value,
          "constant error", NULL },
        { NULL }
};

void __attribute__((__nonnull__))
ls2_add_const_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("const",
                                "Parameters to the constant value error model",
                                "Parameters to the constant value error model",
                                NULL, NULL);
     g_option_group_add_entries(group, const_parameters);
     g_option_context_add_group(context, group);
}


static VECTOR error_v;


void
const_setup(const vector2 *anchors __attribute__((__unused__)),
            size_t nanchors __attribute__((__unused__)))
{
    error_v = VECTOR_BROADCASTF((float) const_arguments.value);
}


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(1,8)))
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
