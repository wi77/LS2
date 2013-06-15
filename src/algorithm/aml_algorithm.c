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

#include "util/util_circle.c"

/*******************************************************************
 ***
 ***   Adapted Multilateration (AML)
 ***
 *******************************************************************/

/* @algorithm_name: Adapted Multilateration (AML) */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
aml_run (const VECTOR* vx, const VECTOR* vy,
         const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy)
{
    VECTOR intersectionx[] = { zero, zero };
    VECTOR intersectiony[] = { zero, zero };
    VECTOR ci, cj;
    VECTOR icount = two;
    
    int is;
    
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        int c = 0;
        float ix[2], iy[2];
        int i,j;
        while (c < 100) {
            i = rand()%(int)num_anchors;
            j = rand()%(int)num_anchors;
            if (j == i) continue;
            is = (int)circle_get_intersection(vx[i][ii],vy[i][ii],vx[j][ii],vy[j][ii], r[i][ii], r[j][ii],ix,iy);
            if (is==2) {
                intersectionx[0][ii] = ix[0];
                intersectionx[1][ii] = ix[1];
                intersectiony[0][ii] = iy[0];
                intersectiony[1][ii] = iy[1];
                ci[ii] = (float)i;
                cj[ii] = (float)j;
                break;
            }
            c++;
        }
    }
     
    /*
    // find two circles (random), which intersect in one or two points
    for (size_t i = 0; i < num_anchors-1; i++) {
        for (size_t j = i+1; j < num_anchors; j++) {
            VECTOR isectx_tmp[2];
            VECTOR isecty_tmp[2];
            VECTOR icount_tmp = vcircle_intersections(vx[i], vy[i], vx[j], vy[j], r[i], r[j], isectx_tmp, isecty_tmp);

            VECTOR mask = VECTOR_AND(VECTOR_EQ(icount, zero), VECTOR_GT(icount_tmp, zero));
            icount = VECTOR_BLENDV(icount, icount_tmp, mask);
            float fi = (float) i; float fj = (float) j;
            ci = VECTOR_BLENDV(ci, VECTOR_BROADCAST(&fi), mask);
            cj = VECTOR_BLENDV(cj, VECTOR_BROADCAST(&fj), mask);
            intersectionx[0] = VECTOR_BLENDV(intersectionx[0], isectx_tmp[0], mask);
            intersectionx[1] = VECTOR_BLENDV(intersectionx[1], isectx_tmp[1], mask);
            intersectiony[0] = VECTOR_BLENDV(intersectiony[0], isecty_tmp[0], mask);
            intersectiony[1] = VECTOR_BLENDV(intersectiony[1], isecty_tmp[1], mask);

            if (VECTOR_TEST_ALL_ONES(VECTOR_GT(icount, zero)))
                goto endloop;
        }
    }
    endloop: ;
    */
    
    
    // Refinement anchors.
    size_t rcount = num_anchors - 2;
    VECTOR rvx[rcount];
    VECTOR rvy[rcount];
    VECTOR rr[rcount];
       
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        // store unused anchors for refinement step
        for (size_t k = 0, n = 0; k < num_anchors; k++) {
            if (k != ci[ii] && k != cj[ii]) {
                rvx[n][ii] = vx[k][ii];
                rvy[n][ii] = vy[k][ii];
                rr[n][ii] = r[k][ii];
                n++;
            }
        }
    }
    
    // Starting point for refinements.
    VECTOR p_intersectx = zero;
    VECTOR p_intersecty = zero;
    
    // if  more than one intersection point, do elimination step
    VECTOR delta_p1 = VECTOR_ABS(distance(intersectionx[0], intersectiony[0], rvx[0], rvy[0]) - rr[0]);
    VECTOR delta_p2 = VECTOR_ABS(distance(intersectionx[1], intersectiony[1], rvx[0], rvy[0]) - rr[0]);
    VECTOR emask = VECTOR_AND(VECTOR_GT(icount, one), VECTOR_GE(delta_p1, delta_p2));
    p_intersectx = VECTOR_BLENDV(intersectionx[0], intersectionx[1], emask);
    p_intersecty = VECTOR_BLENDV(intersectiony[0], intersectiony[1], emask);
    
    // If not intersection, return NAN.
    VECTOR imask = VECTOR_GT(icount, zero);
    p_intersectx = VECTOR_BLENDV(VECTOR_BROADCASTF(NAN), p_intersectx, imask);
    p_intersecty = VECTOR_BLENDV(VECTOR_BROADCASTF(NAN), p_intersecty, imask);
    
    // do first estimation step and refinement
    VECTOR mzero = VECTOR_BROADCASTF(-0.0f);
    for (size_t k = 0; k < rcount; k++) {
        VECTOR pvx = p_intersectx - rvx[k];
        VECTOR pvy = p_intersecty - rvy[k];
        VECTOR dist = VECTOR_SQRT(pvx * pvx + pvy * pvy);
        VECTOR w = rr[k] / dist;
        p_intersectx += (w * pvx + VECTOR_XOR(pvx, mzero)) / two;
        p_intersecty += (w * pvy + VECTOR_XOR(pvy, mzero)) / two;
    }    

    *resx = p_intersectx;
    *resy = p_intersecty;
}
