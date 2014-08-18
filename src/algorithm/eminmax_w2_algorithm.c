/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013  Heiko Will, Marcel Kyas, Thomas Hillebrandt,
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

#include <float.h>
#include <math.h>

#include <glib.h>

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Extended Min Max (W2).
 ***
 *******************************************************************/

/* @algorithm_name: Extended MinMax W2 */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
eminmax_w2_run (const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy) {
    VECTOR east  = vx[0] + r[0];
    VECTOR west  = vx[0] - r[0];
    VECTOR south = vy[0] + r[0];
    VECTOR north = vy[0] - r[0];
    for (size_t i = 1; i < num_anchors; i++) {
        east  = VECTOR_MIN(east,  vx[i] + r[i]);
        west  = VECTOR_MAX(west,  vx[i] - r[i]);
        south = VECTOR_MIN(south, vy[i] + r[i]);
        north = VECTOR_MAX(north, vy[i] - r[i]);
    }
    
    VECTOR corners_x[4];
    VECTOR corners_y[4];
    
    corners_x[0] = west;
    corners_y[0] = north;
    corners_x[1] = east;
    corners_y[1] = north;
    corners_x[2] = west;
    corners_y[2] = south;
    corners_x[3] = east;
    corners_y[3] = south;
        
    VECTOR weights[4];
    for (int j = 0; j < 4; j++) {
        weights[j] = VECTOR_ZERO();
        for (size_t i = 0; i < num_anchors; i++) {		
            weights[j] += (distance(corners_x[j], corners_y[j], vx[i], vy[i]) - r[i]) * (distance(corners_x[j], corners_y[j], vx[i], vy[i]) - r[i]); // W-2 (best for complete scenario)
        }
       	weights[j] = one / VECTOR_MAX(weights[j], VECTOR_BROADCASTF(FLT_EPSILON)); // weight is zero if there are no errors.
    }

    VECTOR ptsx = VECTOR_ZERO(), ptsy = VECTOR_ZERO(), mass = VECTOR_ZERO();
    for (int i = 0; i < 4;  i++) {
        ptsx += corners_x[i] * weights[i];
        ptsy += corners_y[i] * weights[i];
        mass += weights[i];
    }
    *resx = ptsx / mass;
    *resy = ptsy / mass;
}
