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

/* @algorithm_name: CRLB (Zhao Yubin) */

static float crlb_zhao_variance = 25.0f * 25.0f;

#if HAVE_POPT_H
struct poptOption crlb_zhao_arguments[] = {
	{ "crlb-zhao-variance", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &crlb_zhao_variance, 0,
          "Variance of the error model", NULL },
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
crlb_zhao_run(const vector2 *anchor, const size_t num_anchors,
              const vector2 *location)
{
    // Calculate the CRLB
    float numer = 0.0;
    float denom = 0.0;
    for (size_t i = 0; i < num_anchors; i++) {
	const float X_i = location->x - anchor[i].x;
	const float Y_i = location->y - anchor[i].y;
	numer += 1 / crlb_zhao_variance;
	for (size_t j = 0; j < num_anchors; j++) {
	    const float X_j = location->x - anchor[j].x;
	    const float Y_j = location->y - anchor[j].y;
	    denom += (X_i * X_i * Y_j * Y_j - X_i * X_j * Y_i * Y_j) /
                ((X_i * X_i + Y_i * Y_i) * (X_j * X_j + Y_j * Y_j) *
		 (crlb_zhao_variance * crlb_zhao_variance));
	}
    }
    float crlb = numer / denom;
    return crlb;
}
