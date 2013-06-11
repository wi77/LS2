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

/* @algorithm_name: RWGH */

/*******************************************************************
 ***
 ***   RWGH;
 *** 
 *******************************************************************/

#ifndef RWGH_ALGORITHM_C_INCLUDED
#define RWGH_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include "util/util_math.c"
#include "algorithm/nllsq_algorithm.c"

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
rwgh_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    if (width==height){};

    int s = 3;
    int M = (int)no_anchors;
    
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        // calculate all k-permutations for k = s to 
        int int_count = 0;
        int ccount = 0;
        for (int k = s; k <= M; k++) {
            ccount += binom(M, k);
        }

        float intermediatePositions_x[ccount];
        float intermediatePositions_y[ccount];
        float resError[ccount];


        for (int k = s; k <= M; k++) {
            int bino = binom(M, k);
            int permutations[k];
            VECTOR tmpRanges[k];
            VECTOR tmpAnchors_x[k];
            VECTOR tmpAnchors_y[k];
            VECTOR tresx, tresy;
            // initialisation for calculating k-permutations
            for (int i = 0; i < k; i++) {
                permutations[i] = i;
            }

            // build all k-permutations
            for (int i = 0; i < bino; i++) {

                // calculate intermediate position estimate using non-linear
                // least squares multilateration
                for (int h = 0; h < k; h++) {
                    tmpAnchors_x[h][ii] = vx[permutations[h]][ii];
                    tmpAnchors_y[h][ii] = vy[permutations[h]][ii];
                    tmpRanges[h][ii] = r[permutations[h]][ii];
                }

                nllsq_run(tmpAnchors_x, tmpAnchors_y, tmpRanges, (size_t)k, width, height, &tresx, &tresy);
                resError[int_count] = calculate_residual_error((size_t)k, tmpAnchors_x, tmpAnchors_y, tmpRanges, tresx,tresy)[ii] / (float)k;
                intermediatePositions_x[int_count] = tresx[ii];
                intermediatePositions_y[int_count] = tresy[ii];
                int_count++;
                // build next permutation
                if (i == bino - 1) {
                    break;
                }
                int j = k - 1;
                while (j >= 0) {
                    if (!incCounter(permutations, j, M, k)) {
                        break;
                    }
                    j--;
                }
                for (int l = j+1; l < k; l++) {
                    permutations[l] = permutations[l-1] + 1;
                }
            }
        }

        float res = 0;
        float x_x = 0.0f;
        float x_y = 0.0f;

        for (int i = 0; i < int_count; i++) {
            x_x += intermediatePositions_x[i] * (1/resError[i]);
            x_y += intermediatePositions_y[i] * (1/resError[i]);
            res += 1/resError[i];
        }
        x_x /= res;
        x_y /= res;
        (*resx)[ii] = x_x;
        (*resy)[ii] = x_y;
    }
}

#endif
