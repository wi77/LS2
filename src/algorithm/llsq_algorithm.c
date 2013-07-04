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
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Linear Least Squares
 ***
 *******************************************************************/
 
 /* @algorithm_name: Linear Least Squares */

#ifndef ALGORITHM_LLSQ_C_INCLUDED
#define ALGORITHM_LLSQ_C_INCLUDED 1

#include "util/util_matrix.c"
#include "util/util_triangle.c"
 
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
llsq_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
         size_t num_anchors,
         int width __attribute__((__unused__)),
         int height __attribute__((__unused__)),
         VECTOR *restrict resx,
         VECTOR *restrict resy)
{
    // Solve equation of form W*A*x = W*b
    assert(num_anchors > 0);
    const size_t m = num_anchors-1U;
    VECTOR b[m*1];
    VECTOR a[m*2];
    VECTOR mAT[2*m];
    VECTOR tmp[2*m];
    VECTOR tmp2[2*2];
    VECTOR tmp3[2*2];
    VECTOR tmp4[2*1];    

    for (size_t i = 0; i < m; i++) {
        a[i*2+0] = (vx[i] - vx[m]);
        a[i*2+1] = (vy[i] - vy[m]);
        b[i*1+0] = half * (vx[i] * vx[i] - vx[m] * vx[m] +
                    vy[i] * vy[i] - vy[m] * vy[m] +
                    r[m] * r[m] - r[i] * r[i]);
    }
    
    // Solve with closed form solution: x = (A^T*A)^-1 * A^T * b
        transpose(a,mAT,m,2);       
        times(mAT,2,m,a,2,tmp2);
        invert_2x2(tmp2,tmp3);
        times(tmp3,2,2,mAT,m,tmp);        
        times(tmp,2,m,b,1,tmp4);
        
        *resx = tmp4[0];
        *resy = tmp4[1];
}
#endif
