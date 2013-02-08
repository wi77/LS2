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

/* @algorithm_name: CRLB (Qi & Kobayashi) */

static float crlb_malaney_sigma = 3.0f;

static float crlb_malaney_pathloss = 2.0f;

static float crlb_malaney_scale = 1.0f;

#if HAVE_POPT_H
struct poptOption crlb_malaney_arguments[] = {
        { "crlb-malaney-noise", 'S',
          POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &crlb_malaney_sigma, 0,
          "noise variable", NULL },
	{ "crlb-malaney-exponent", 'n',
          POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &crlb_malaney_pathloss, 0,
          "path loss exponent", NULL },
	{ "crlb-malaney-scale", 'Z',
          POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &crlb_malaney_scale, 0,
          "scale factor for root mean square error", NULL },
        POPT_TABLEEND
};
#endif

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
        (10 * crlb_malaney_pathloss) / (crlb_malaney_sigma * logf(10.0f));

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
