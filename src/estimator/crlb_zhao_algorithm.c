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

#ifdef HAVE_CONFIG_H
#  include <ls2/ls2-config.h>
#endif

#include <glib.h>

#include "crlb_zhao_algorithm.h"

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Cramer Rao Lower Bound by Qi and Kobayashi
 ***
 *******************************************************************/

/* @algorithm_name: CRLB (Zhao Yubin) */

typedef struct crlb_zhao_arguments {
    double variance;
} crlb_zhao_arguments;


crlb_zhao_arguments ls2_crlb_zhao_arguments;

void __attribute__((__nonnull__))
ls2_crlb_zhao_init_arguments(struct crlb_zhao_arguments *arguments)
{
	arguments->variance = 25.0 * 25.0;
}


GOptionEntry crlb_zhao_parameters[] = {
	{ "crlb-zhao-variance", 0, 0, G_OPTION_ARG_DOUBLE,
          &ls2_crlb_zhao_arguments.variance,
          "Variance of the error model", NULL },
        { NULL }
};



void __attribute__((__nonnull__))
ls2_add_crlb_zhao_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("crlb-zhao",
                                "Parameters to the CRLB model of Zhao",
                                "Parameters to the CRLB model of Zhao",
                                NULL, NULL);
     g_option_group_add_entries(group, crlb_zhao_parameters);
     g_option_context_add_group(context, group);
}



/*!
 * Do not use an error model with this one, it computes nonesense in this
 * situation! Use a constant error of 0.
 *
 * This algorithm calculates the lower bound of the mean-squared
 * estimation error.
 */

#include <float.h>

static inline float __attribute__((__always_inline__,__pure__))
crlb_zhao_run(const vector2 *anchor, const size_t num_anchors,
              const vector2 *location)
{
    const float variance = (float) ls2_crlb_zhao_arguments.variance;
    // Calculate the CRLB
    float numer = 0.0f;
    float denom = 0.0f;
    for (size_t i = 0; i < num_anchors; i++) {
	const float X_i = location->x - anchor[i].x;
	const float Y_i = location->y - anchor[i].y;
	numer += 1.0f / variance;
	for (size_t j = 0; j < num_anchors; j++) {
	    const float X_j = location->x - anchor[j].x;
	    const float Y_j = location->y - anchor[j].y;
	    denom += (X_i * X_i * Y_j * Y_j - X_i * X_j * Y_i * Y_j) /
                ((X_i * X_i + Y_i * Y_i) * (X_j * X_j + Y_j * Y_j) *
		 (variance * variance));
	}
    }
    float crlb = numer / denom;
    return crlb;
}
