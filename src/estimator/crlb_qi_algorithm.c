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

#include "crlb_qi_algorithm.h"

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

/* @algorithm_name: CRLB (Qi, Kobayashi, Suda) */



// That is 5 million divided by square root of 3, as found in Qi's
// thesis (p. 18, eq. 2.29)
static double crlb_qi_beta = 2886751.345948129;

static double crlb_qi_su = 0.1;

static double crlb_qi_scale = 1.0;

GOptionEntry crlb_qi_arguments[] = {
        { "bandwidth", 0, 0, G_OPTION_ARG_DOUBLE,
          &crlb_qi_beta,
          "effective bandwidth of the signal waveform", NULL },
	{ "scale", 0, 0, G_OPTION_ARG_DOUBLE,
          &crlb_qi_scale,
          "scale factor for root mean square error", NULL },
	{ "unit", 0, 0, G_OPTION_ARG_DOUBLE,
          &crlb_qi_su,
          "length of a simulation unit in meters", NULL },
        { NULL }
};

void __attribute__((__nonnull__))
ls2_add_crlb_qi_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("crlb-qi",
                                "Parameters to the CRLB model of Qi and Kobayashi",
                                "Parameters to the CRLB model of Qi and Kobayashi",
                                NULL, NULL);
     g_option_group_add_entries(group, crlb_qi_arguments);
     g_option_context_add_group(context, group);
}


#include "algorithm/trilateration_algorithm.c"

/*!
 * 
 *
 * Do not use an error model with this one, it computes nonesense in this
 * situation! Use a constant error of 0.
 *
 * This algorithm calculates the lower bound of the mean-squared
 * estimation error.
 *
 * See Y. Qi, H. Kobayashi, H. Suda. Analysis of Wireless Geolocation
 * in a Non-Line-of-Sight Environment. IEEE Transactions of Wireless
 * Communications, Vol. 5, No. 3, March 2006, pp. 672-681.
 *
 * This implements Equation (35) of the above paper. We assume that
 * all nodes have line of sight.
 */
static inline float __attribute__((__always_inline__,__pure__))
crlb_qi_run(const vector2 *anchor, const size_t count, const vector2 *location)
{
    const float pi = (float) M_PI;
    // speed of light in su / s, assumes su = 1dm
    const float c = (float) (299792458.0 / crlb_qi_su);
    // Constant alpha of the paper
    const float alpha =
         (float) ((c * c) / (8.0 * pi * pi * crlb_qi_beta * crlb_qi_beta));

    // Calculate the CRLB
    float numer = 0.0f;
    float denom = 0.0f;
    for (size_t i = 0; i < count; i++) {
	numer += distance_v(&(anchor[i]), location);
	for (size_t j = 0; j < i; j++) {
	    float phi_i = atanf((location->y - anchor[i].y) / (location->x - anchor[i].x));
	    float phi_j = atanf((location->y - anchor[j].y) / (location->x - anchor[j].x));
	    float a = sinf(phi_i - phi_j);
	    denom += distance_v(location, &(anchor[i])) *
                       distance_v(location, &(anchor[j])) * a * a;
	}
    }
    return alpha * numer / denom;
}
