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

/********************************************************************
 **
 **  This file is made only for including in the LS² project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/* @algorithm_name: BILATERATION */

/*******************************************************************
 ***
 ***   BILATERATION;
 *** 
 *******************************************************************/

#ifndef BILATERATION_ALGORITHM_C_INCLUDED
#define BILATERATION_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include "util/util_math.c"

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
bilateration_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    if (width==height){};
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        // step 1: calculate circle intersections with approximation
        int n = (int)no_anchors;
        int bino = binom(n, 2);
        float intersections_x[bino][2];
        float intersections_y[bino][2];
        size_t intersection_count[bino];

        // build all k-permutations
        int num = 0;
        size_t is = 0;
        for (int i = 0; i < n-1; i++) {
            for (int j = i+1; j < n; j++) {
                is = circle_get_intersection (vx[i][ii], vy[i][ii], vx[j][ii], vy[j][ii], r[i][ii], r[j][ii], &(intersections_x[num][0]), &(intersections_y[num][0]));
                intersection_count[num] = is;
                // no intersection => try to get approximated intersection
                if (!is) {
                    float pft_x,pft_y;
                    is = circle_get_approx_intersection2(vx[i][ii], vy[i][ii], vx[j][ii], vy[j][ii], r[i][ii], r[j][ii], &pft_x, &pft_y);
                    if (is) {
                        intersections_x[num][0] = pft_x;
                        intersections_y[num][0] = pft_y;
                        intersection_count[num] = 1;
                    }
                }
                if (is) {
                    num++;
                }
            }
        }

        float listT_x[num*2];
        float listT_y[num*2];
        int lcount = 0;
        
        // List of points finally used for localization
        for (int i = 0; i < num; i++) {
           // No decision on single intersections
           // Point2d[] currentI = intersections[i];
            if (intersection_count[i] == 1) {
                listT_x[lcount] = intersections_x[i][0];
                listT_y[lcount] = intersections_y[i][0];
                lcount++;
                continue;
            }
            double psi = 0;
            double phi = 0;
            for (int j = 0; j < num; j++) {
                if (i != j) {
                    double deltaPsi;
                    double deltaPhi;
                    //Point2d[] compareI = intersections[j];
                    if (intersection_count[j] == 1) {
                        deltaPsi = distance_squared_s(intersections_x[i][0], intersections_y[i][0], intersections_x[j][0], intersections_y[j][0]);
                        deltaPhi = distance_squared_s(intersections_x[i][1], intersections_y[i][1], intersections_x[j][0], intersections_y[j][0]);
                    } else {
                        deltaPsi =  fminf(distance_squared_s(intersections_x[i][0], intersections_y[i][0], intersections_x[j][0], intersections_y[j][0]),
                                        distance_squared_s(intersections_x[i][0], intersections_y[i][0], intersections_x[j][1], intersections_y[j][1]));
                        deltaPhi =  fminf(distance_squared_s(intersections_x[i][1], intersections_y[i][1], intersections_x[j][0], intersections_y[j][0]),
                                        distance_squared_s(intersections_x[i][1], intersections_y[i][1], intersections_x[j][1], intersections_y[j][1]));
                    }
                    psi += deltaPsi;
                    phi += deltaPhi;
                }
            }
            if (psi < phi) {
                listT_x[lcount] = intersections_x[i][0];
                listT_y[lcount] = intersections_y[i][0];
                lcount++;
            } else {
                listT_x[lcount] = intersections_x[i][1];
                listT_y[lcount] = intersections_y[i][1];
                lcount++;
            }
        }
        
        // Final position is centroid of all left points
        float mass[lcount];
        for (int jj=0;jj<lcount;jj++)
            mass[jj]=1.0f;
        center_of_mass(lcount, listT_x, listT_y, mass, &((*resx)[ii]), &((*resy)[ii]));
    }
}

#endif
