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

/* @algorithm_name: RLSM */

/*******************************************************************
 ***
 ***   RLSM;
 *** 
 *******************************************************************/

#ifndef RLSM_ALGORITHM_C_INCLUDED
#define RLSM_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include "util/util_math.c"
#include "util/util_points.c"
#include "algorithm/nllsq_algorithm.c"

static inline int
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
compare_floats2 (const void *a_, const void *b_) 
{
  const float *a = a_;
  const float *b = b_;

  return *a < *b ? -1 : *a > *b;
}

static inline int 
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
toIndex(int row, int col, int N) {
        int idx = -1;
        if (row < col) {
            idx = row * (N-1) - (row-1) * ((row-1) + 1)/2 + col - row - 1;
        } else if (col < row) {
            idx = col * (N-1) - (col-1) * ((col-1) + 1)/2 + row - col - 1;
        }
        return idx;
    }
    
static inline int 
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
robust_filter(int const count ,float *restrict resx,float *restrict resy) {
        int N = binom(count, 2);
        float v[N];
        float vC[N];
        for (int i = 0; i < count-1; i++) {
            for (int j = i+1; j < count; j++) {
                int idx = toIndex(i, j, count);
                v[idx] = distance_s(resx[i],resy[i],resx[j],resy[j]);
                vC[idx] = v[idx];
            }
        }

        qsort(vC,(size_t)N,sizeof(float),compare_floats2);
        float MEDV = 2 * vC[count/2+1];
        float filtered_x[count];
        float filtered_y[count];
        int filtered_count = 0;
        for (int i = 0; i < count; i++) {
            int dropCounter = 0;
            for (int j = 0; j < count; j++) {
                if (i == j) {
                    continue;
                }
                double dist = v[toIndex(i, j, count)];
                if (dist >= MEDV) {
                    dropCounter++;
                }
            }
            if (dropCounter <= count/2) {
                filtered_x[filtered_count] = resx[i];
                filtered_y[filtered_count] = resy[i];
                filtered_count++;
            }
        }
        memcpy(resx,filtered_x,sizeof(float)*(unsigned int)filtered_count);
        memcpy(resy,filtered_y,sizeof(float)*(unsigned int)filtered_count);
        return filtered_count;
    }    


static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
rlsm_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
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
                    tmpAnchors_x[h][i%VECTOR_OPS] = vx[permutations[h]][ii];
                    tmpAnchors_y[h][i%VECTOR_OPS] = vy[permutations[h]][ii];
                    tmpRanges[h][i%VECTOR_OPS] = r[permutations[h]][ii];
                }

                if (((i+1)%VECTOR_OPS)==0 || i == bino - 1){
                    llsq_run(tmpAnchors_x, tmpAnchors_y, tmpRanges, (size_t)k, width, height, &tresx, &tresy);
                    for (int jj=0; jj <= i%VECTOR_OPS; jj++){
                        if (!isnan(tresy[jj])){
                            intermediatePositions_x[int_count] = tresx[jj];
                            intermediatePositions_y[int_count] = tresy[jj];
                            int_count++;
                        }
                    }
                }                

                
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

        float x_x = 0.0f;
        float x_y = 0.0f;

        // apply robust median filter to this array of intermediate
        // position estimates
        if (int_count > 1) {
            int_count = robust_filter(int_count, intermediatePositions_x, intermediatePositions_y);
        }
        
        float weights[int_count];
        for (int jj=0;jj<int_count;jj++)
            weights[jj]=1.0f;

        // return geometric median as result
        point_geometric_median(int_count, intermediatePositions_x,
                       intermediatePositions_y,
                       weights,
                       &x_x, &x_y);
      
        (*resx)[ii] = x_x;
        (*resy)[ii] = x_y;
    }
}

#endif
