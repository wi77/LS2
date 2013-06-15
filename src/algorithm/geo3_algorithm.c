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

#ifndef GEO3_ALGORITHM_C_INCLUDED
#define GEO3_ALGORITHM_C_INCLUDED 1

#include "util/util_median.c"
#include "util/util_points_opt.c"

/*******************************************************************
 ***
 ***   Geolateration
 ***
 *******************************************************************/

/* @algorithm_name: Geolateration (3 anchors only) */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
geo3_run(const const VECTOR* vx, const VECTOR* vy,
         const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy)
{
    if (num_anchors<3) return;
    
    // Step 0: Precalculate mandatory values
    
    // real circle intersections
    VECTOR _rc[3];
    VECTOR _ix[6];
    VECTOR _iy[6];
    _rc[0] = vcircle_intersections(vx[0], vy[0], vx[1], vy[1], r[0], r[1], _ix, _iy);
    _rc[1] = vcircle_intersections(vx[0], vy[0], vx[2], vy[2], r[0], r[2], &_ix[2], &_iy[2]);
    _rc[2] = vcircle_intersections(vx[1], vy[1], vx[2], vy[2], r[1], r[2], &_ix[4], &_iy[4]);
    
    // anchor center
    VECTOR icr = triangle_icr(vx[0], vy[0], vx[1], vy[1], vx[2], vy[2]);
    
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        size_t rc, sum = 0;
        float ix[8];    // REMARK: Indices 6 & 7 are not used.
        float iy[8];

        // step 1: calculate approximated circle intersections where necessary
        rc = (size_t) _rc[0][ii]; // is a non-negative integer.
        ix[0] = _ix[0][ii]; ix[1] = _ix[1][ii];
        iy[0] = _iy[0][ii]; iy[1] = _iy[1][ii];
        if (rc == 0) {
            // no intersection => try to get approximated intersection
            rc = circle_get_approx_intersection(vx[0][ii], vy[0][ii], vx[1][ii], vy[1][ii], r[0][ii],  r[1][ii], ix, iy);
            if (rc == 0) {
                 (*resx)[ii] = NAN;
                 (*resy)[ii] = NAN;
                 continue;           
            };
        }
        sum += rc;
        
        rc = (size_t) _rc[1][ii]; // is a non-negative integer.
        ix[sum] = _ix[2][ii]; ix[sum+1] = _ix[3][ii];
        iy[sum] = _iy[2][ii]; iy[sum+1] = _iy[3][ii];
        if (rc == 0) {
            // no intersection => try to get approximated intersection
            rc = circle_get_approx_intersection(vx[0][ii], vy[0][ii], vx[2][ii], vy[2][ii], r[0][ii],  r[2][ii], &ix[sum], &iy[sum]);
            if (rc == 0) {
                 (*resx)[ii] = NAN; 
                 (*resy)[ii] = NAN;
                 continue;           
            };
        }
        sum += rc;
        
        rc = (size_t) _rc[2][ii]; // is a non-negative integer.
        ix[sum] = _ix[4][ii]; ix[sum+1] = _ix[5][ii];
        iy[sum] = _iy[4][ii]; iy[sum+1] = _iy[5][ii];
        if (rc == 0) {
            // no intersection => try to get approximated intersection
            rc = circle_get_approx_intersection(vx[1][ii], vy[1][ii], vx[2][ii], vy[2][ii], r[1][ii],  r[2][ii], &ix[sum], &iy[sum]);
            if (rc == 0) {
                 (*resx)[ii] = NAN;
                 (*resy)[ii] = NAN;
                 continue;           
            };
        }
        sum +=rc;
        
        // Precalculate distances
#ifdef __AVX__
        union {
            VECTOR v;
            float f[8];
        } distances[sum];
        VECTOR xtargets = _mm256_load_ps(ix);
        VECTOR ytargets = _mm256_load_ps(iy);
        for(size_t i = 0; i < sum; i++) {
            VECTOR px = VECTOR_BROADCAST(&ix[i]);
            VECTOR py = VECTOR_BROADCAST(&iy[i]);
            distances[i].v = distance(px, py, xtargets, ytargets);
            for(size_t j = sum; j < 8; j++)
                distances[i].f[j] = FLT_MAX;
        }
#else
        union {
            VECTOR v[2];
            float f[8];
        } distances[6];
        VECTOR xtargets_lower = _mm_load_ps(ix);
        VECTOR ytargets_lower = _mm_load_ps(iy);
        VECTOR xtargets_upper = _mm_load_ps(&ix[4]);
        VECTOR ytargets_upper = _mm_load_ps(&iy[4]);
        xtargets_upper = _mm_movelh_ps(xtargets_upper, xtargets_upper);
        ytargets_upper = _mm_movelh_ps(ytargets_upper, ytargets_upper);
        
        for(size_t i = 0; i < sum; i+=2) {
            VECTOR px = VECTOR_BROADCAST(&ix[i]);
            VECTOR py = VECTOR_BROADCAST(&iy[i]);
            VECTOR px2 = VECTOR_BROADCAST(&ix[i + 1]);
            VECTOR py2 = VECTOR_BROADCAST(&iy[i + 1]);
            
            distances[i].v[0] = distance(px, py, xtargets_lower, ytargets_lower);
            distances[i + 1].v[0] = distance(px2, py2, xtargets_lower, ytargets_lower);
                        
            px = _mm_movelh_ps(px, px2);
            py = _mm_movelh_ps(py, py2);
            distances[i].v[1] = distance(px, py, xtargets_upper, ytargets_upper);
            
            distances[i + 1].f[4] = distances[i].f[6];
            distances[i + 1].f[5] = distances[i].f[7];
            
            for(size_t j = sum; j < 8; j++) {
                distances[i].f[j] = FLT_MAX;
                distances[i + 1].f[j] = FLT_MAX;
            }
        }
#endif
            
        // step 2: if there are 3 points which are very close together
        //          => no ranging error, take one of them as result
        int found = 0;
        VECTOR max_distance = VECTOR_BROADCASTF(0.1f);
        for (size_t i = 0; i < sum; i++) {     
#ifdef __AVX__         
            VECTOR tmp = VECTOR_AND(one, VECTOR_LT(distances[i].v, max_distance));
            tmp = _mm256_hadd_ps(tmp, tmp);
            tmp = _mm256_hadd_ps(tmp, tmp);
            tmp = _mm256_hadd_ps(tmp, tmp);
#else
            VECTOR tmp = VECTOR_AND(one, VECTOR_LT(distances[i].v[0], max_distance));
            tmp += VECTOR_AND(one, VECTOR_LT(distances[i].v[1], max_distance));
            tmp = _mm_hadd_ps(tmp, tmp);
            tmp = _mm_hadd_ps(tmp, tmp);
#endif
            if(tmp[0] >= 3) {
                (*resx)[ii] = ix[i];
                (*resy)[ii] = iy[i];
                found = 1;
                break;
            }
        }
        if(found)
            continue;

        // step 3: calculate triangle with minimum perimeter, also check for
        //         minimum triangle which lies inside of all circles.
        int min1 = -1, min2 = -1, min3 = -1;
        int minIn1 = -1, minIn2 = -1, minIn3 = -1;
        float minV = FLT_MAX, minInV = FLT_MAX;
        float anchorsx[] = {vx[0][ii], vx[1][ii], vx[2][ii], 0};
        float anchorsy[] = {vy[0][ii], vy[1][ii], vy[2][ii], 0};

        // Precalculate point in circles booleans.
        __m128 vanchorsx = _mm_load_ps(anchorsx);
        __m128 vanchorsy = _mm_load_ps(anchorsy);
        __m128 vmaxdists = _mm_setr_ps(r[0][ii] + 0.01f, r[1][ii] + 0.01f, r[2][ii] + 0.01f, FLT_MAX);
        int pic[sum];
        for(size_t i = 0; i < sum; i++)
            pic[i] = circle_point_in_circles_opt128(vanchorsx, vanchorsy, vmaxdists, ix[i], iy[i]);
        
        
        for (size_t i = 0; i < sum - 2; i++) {
            for (size_t j = i + 1; j < sum - 1; j++) {
                for (size_t k = j + 1; k < sum; k++) {
                    float tmp = distances[i].f[j] + distances[i].f[k] + distances[j].f[k];
                    // test if new minimum triangle found
                    if (tmp < minV) {
                        minV = tmp;
                        min1 = (int) i;
                        min2 = (int) j;
                        min3 = (int) k;
                    }
                    
                    // test for minimum triangle which lies in all circles
                    if (pic[i] + pic[j] + pic[k] == 3) {
                        if (tmp < minInV) {
                            minInV = tmp;
                            minIn1 = (int) i;
                            minIn2 = (int) j;
                            minIn3 = (int) k;
                        }
                    }
                }
            }
        }


        // REMARK: the first implementation of this algorithm stopped at this
        //         point and returned the center of gravity of the minimum
        //         triangle as result (without finding minimum triangle which
        //         lies in all circles).

        // step 4: (small optimization) Test if center of gravity of the min
        //         triangle lies in the center of the triangle formed by the
        //         anchor nodes. Therefore, calculate incircle radius of the
        //         triangle reduced by some factor.
        
        // Calculate center of mass here (not in util_points) to optimize performance.
        float pcomx =  (ix[min1] + ix[min2] + ix[min3]) / 3;
        float pcomy = (iy[min1] + iy[min2] + iy[min3]) / 3;
        float panchorcomx = (anchorsx[0] + anchorsx[1] + anchorsx[2]) / 3;
        float panchorcomy = (anchorsy[0] + anchorsy[1] + anchorsy[2]) / 3;

        if (distance_s(panchorcomx, panchorcomy,pcomx,pcomy) <= icr[ii]) {
            (*resx)[ii] = pcomx;
            (*resy)[ii] = pcomy;
            continue;
        }

        // step 5: test if there is a minimum triangle which lies in all
        //         circles. If it's not the "real" minimum triangle then
        //         do some further optimization.
        if (minIn1 != -1) {
            int same = (min1 == minIn1 && min2 == minIn2 && min3 == minIn3);
            if (!same) {
                // calculate the area of both tiangles
                float ta1 = triangle_area_s(ix[min1], iy[min1], ix[min2], iy[min2], ix[min3], iy[min3]);
                float ta2 = triangle_area_s(ix[minIn1], iy[minIn1], ix[minIn2], iy[minIn2], ix[minIn3], iy[minIn3]);
                
                // if the area of the second triangle is roughly equal to the
                // area of the minimum triangle, use the second triangle.
                if (ta2 < 2.5 * ta1) {
                    min1 = minIn1;
                    min2 = minIn2;
                    min3 = minIn3;
                } else {
                    // build weights with triangle areas
                    ta1 = 1.0f / (ta1);
                    ta2 = 1.0f / (ta2);
                    
                    // calculate weighted geometric median as result
                    float ccmx[] = {ix[min1], ix[min2], ix[min3], ix[minIn1], ix[minIn2], ix[minIn3]};
                    float ccmy[] = {iy[min1], iy[min2], iy[min3], iy[minIn1], iy[minIn2], iy[minIn3]};
                    float taw[] = {ta1, ta1, ta1, ta2, ta2, ta2};
                    point_geometric_median(6, ccmx, ccmy, taw, &((*resx)[ii]), &((*resy)[ii]));
                    continue;
                }
            }
        }
        
        __m128 ccmx = _mm_setr_ps(ix[min1], ix[min2], ix[min3], 0.0);
        __m128 ccmy = _mm_setr_ps(iy[min1], iy[min2], iy[min3], 0.0);
        // step 6: calculate geometric median of final triangle as result
        __m128 mass = _mm_setr_ps(1.0, 1.0, 1.0, 0.0);
        point_geometric_median_opt(3, ccmx, ccmy, mass, &((*resx)[ii]), &((*resy)[ii]) );
    }
}

#endif
