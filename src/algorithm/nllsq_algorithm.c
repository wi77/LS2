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

/* @algorithm_name: Nonlinear Least Squares */

/*******************************************************************
 ***
 ***   Nonlinear Least Squares
 ***
 *******************************************************************/


#ifndef NLLSQ_ALGORITHM_C_INCLUDED
#define NLLSQ_ALGORITHM_C_INCLUDED 1

#include "llsq_algorithm.c"

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
nllsq_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors,
          int width, int height, VECTOR *restrict resx, VECTOR *restrict resy)
{
        float load = -1.0F;
        *resx = VECTOR_BROADCAST(&load);
        *resy = VECTOR_BROADCAST(&load);
        // Starting point of optimization
        VECTOR s0x;
        VECTOR s0y;

        // 1a. Set starting point of optimization to linear least squares result
        llsq_run (vx, vy, r, num_anchors, width, height, &s0x, &s0y);
        
        VECTOR e0, e1;
        int iterations = 0;
        float epsilon = 0.000001F;
        VECTOR r0x;
        VECTOR r0y;

        do {
            // initial squared error
            e0 = calculate_residual_error(num_anchors, vx,vy,r,s0x,s0y);

            // 2. Solve equation of form A*x = b, A is A.length x 2
            VECTOR b[num_anchors*1];
            VECTOR a[num_anchors*2];
            for (size_t i = 0; i < num_anchors; i++) {
                VECTOR dist = distance(s0x,s0y,vx[i],vy[i]);
                a[i*2+0] = (s0x - vx[i]) / dist;
                a[i*2+1] = (s0y - vy[i]) / dist;
                b[i*1+0] = (r[i] - dist) + (a[i*2+0] * s0x + a[i*2+1] * s0y);
            }

            // Solve with closed form solution: x = (A^T*A)^-1 * A^T * b
            VECTOR mAT[2*num_anchors];
            VECTOR tmp[2*num_anchors];
            VECTOR tmp2[2*2];
            VECTOR tmp3[2*2];
            VECTOR tmp4[2*1];
     
            transpose(a,mAT,num_anchors,2);       
            times(mAT,2,num_anchors,a,2,tmp2);
            invert_2x2(tmp2,tmp3);
            times(tmp3,2,2,mAT,num_anchors,tmp);        
            times(tmp,2,num_anchors,b,1,tmp4);

            r0x = tmp4[0];
            r0y = tmp4[1];

            // new squared error
            e1 = calculate_residual_error(num_anchors, vx,vy,r,r0x,r0y);
            VECTOR temp = e0 - e1;

            int br=1;
            for (int i=0; i < VECTOR_OPS; i++){
                if ((*resx)[i] == -1.0f) {
                    br = 0;
                    if (temp[i] < epsilon){
                        (*resx)[i] = r0x[i];
                        (*resy)[i] = r0y[i];
                    }        
                }
            }
            if (br) break;
            
            // Set refined position for next step
            s0x = r0x;
            s0y = r0y;
            iterations++;

        } while (iterations <= 100);
        for (int i=0; i < VECTOR_OPS; i++){
            if ((*resx)[i] == -1.0F) {
                 (*resx)[i] = s0x[i];
                 (*resy)[i] = s0y[i];      
            }           
        }
}

#endif
