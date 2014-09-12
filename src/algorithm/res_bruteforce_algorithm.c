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
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Min-Max aka Bounding Box algorithm.
 ***
 *******************************************************************/

/* @algorithm_name: Residual-Bruteforce*/

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
res_bruteforce_run (const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors, size_t width __attribute__((__unused__)), size_t height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy) {
    float min_res;
    int minres_x = 0.0,minres_y = 0.0;
   
    int x,y,a,i;
    for (i=0;i<VECTOR_OPS;i++) {
        min_res = FLT_MAX;
        for (y=0;y<250;y++){
            for (x=0;x<250;x++){
                float res = 0;
                float d = 0;
                for (a=0; a < (int)num_anchors; a++) {
                    d =  sqrtf(((float)x-vx[a][i])*((float)x-vx[a][i])+((float)y-vy[a][i])*((float)y-vy[a][i])) - r[a][i];
                    res += d*d;
                } 
                if (res < min_res) {
                    min_res = res;
                    minres_x = x;
                    minres_y = y;
                }
            }
        }
        (*resx)[i] = (float)minres_x;
        (*resy)[i] = (float)minres_y;   
    }
}
