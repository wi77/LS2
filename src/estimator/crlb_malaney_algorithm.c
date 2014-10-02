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

#include <glib.h>

#include "crlb_malaney_algorithm.h"


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

/* @algorithm_name: CRLB (Malaney) */


typedef struct crlb_malaney_arguments {
    double sigma;
    double pathloss;
    double scale;
} crlb_malaney_arguments;

crlb_malaney_arguments ls2_crlb_malaney_arguments;

void __attribute__((__nonnull__))
ls2_crlb_malaney_init_arguments(struct crlb_malaney_arguments *arguments)
{
     arguments->sigma = 3.0;
     arguments->pathloss = 2.0;
     arguments->scale = 1.0;
}

GOptionEntry crlb_malaney_parameters[] = {
        { "crlb-malaney-noise", 0, 0, G_OPTION_ARG_DOUBLE,
          &(ls2_crlb_malaney_arguments.sigma),
          "noise variable", NULL },
	{ "crlb-malaney-exponent", 0, 0, G_OPTION_ARG_DOUBLE,
          &(ls2_crlb_malaney_arguments.pathloss),
          "path loss exponent", NULL },
	{ "crlb-malaney-scale", 0, 0, G_OPTION_ARG_DOUBLE,
          &(ls2_crlb_malaney_arguments.scale),
          "scale factor for root mean square error", NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_crlb_malaney_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("crlb-malaney",
                                "Parameters to the CRLB model of Malaney",
                                "Parameters to the CRLB model of Malaney",
                                NULL, NULL);
     g_option_group_add_entries(group, crlb_malaney_parameters);
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
crlb_malaney_run (const vector2 *anchor, const size_t num_anchors,
	          const vector2 *location)
{
    const float alpha =
        (float) ((10.0 * ls2_crlb_malaney_arguments.pathloss) / (ls2_crlb_malaney_arguments.sigma * log(10.0)));

    // Calculate the CRLB
    float numer = 0.0f;
    float denom = 0.0f;
    for (size_t i = 0; i < num_anchors - 1; i++) {
        float d = distance_v(&(anchor[i]), location);
	numer += 1/(d * d);
	for (size_t j = i + 1; j < num_anchors; j++) {
            float e = distance_v(&(anchor[j]), location);
	    float phi_i = atanf((location->x - anchor[j].x) / (location->x - anchor[i].x));
	    float phi_j = atanf((location->y - anchor[j].y) / (location->x - anchor[i].y));
	    float a = sinf(phi_i - phi_j);
	    denom += a * a / (d * d * e * e);
	}
    }
    return numer / (alpha * alpha * denom);
}
