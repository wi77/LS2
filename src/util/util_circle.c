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

#ifndef INCLUDED_UTIL_CIRCLE_C
#define INCLUDED_UTIL_CIRCLE_C 1

#include "util/util_vector.c"

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/
 
 /*******************************************************************
 ***
 *** Circle Operations
 ***
 *******************************************************************/

/**
 * Returns the intersections of two circles.
 * @ret: the number of intersection [0,1,2]
 */
static inline size_t
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
circle_get_intersection (const float  p1x, const float  p1y, const float  p2x,
                         const float  p2y, const float r1, const float r2,
                         float *restrict retx, float *restrict rety)
{
    float d = distance_s(p1x,p1y,p2x,p2y);

    // no solutions, the circles are separate || the circles are coincident || no solutions because one circle is contained within the other
    // => infinite number of solutions possible
    if (r1+r2 < d || fabs(r1-r2) > d || d == 0) {
        return 0;
    }

    float a = (r1*r1 - r2*r2 + d*d) / (2.0f * d);
    float v = r1*r1 - a*a;
    float h = sqrtf(v);

    float dx = (p2x - p1x) / d;
    float dy = (p2y - p1y) / d;
    float p3x = p1x + a * dx;
    float p3y = p1y + a * dy;

    dx *= h;
    dy *= h;
    float p4x = p3x + dy;
    float p4y = p3y - dx;

    size_t count = (p4x == p3x && p4y == p3y) ? 1 : 2;

    retx[0] = p4x;
    rety[0] = p4y;
    if (count == 2) {
        p4x = p3x - dy;
        p4y = p3y + dx;
        retx[1] = p4x;
        rety[1] = p4y;
    }
    return count;
}

static inline size_t
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
circle_get_intersectionf (const float  p1x, const float  p1y, const float  p2x,
                         const float  p2y, const float r1, const float r2,
                         double *restrict retx, double *restrict rety)
{
    double d = distance_sf(p1x,p1y,p2x,p2y);

    // no solutions, the circles are separate || the circles are coincident || no solutions because one circle is contained within the other
    // => infinite number of solutions possible
    if (r1+r2 < d || fabs(r1-r2) > d || d == 0) {
        return 0;
    }

    double a = (r1*r1 - r2*r2 + d*d) / (2.0f * d);
    double v = r1*r1 - a*a;
    double h = sqrt(v);

    double dx = (p2x - p1x) / d;
    double dy = (p2y - p1y) / d;
    double p3x = p1x + a * dx;
    double p3y = p1y + a * dy;

    dx *= h;
    dy *= h;
    double p4x = p3x + dy;
    double p4y = p3y - dx;

    size_t count = (p4x == p3x && p4y == p3y) ? 1 : 2;

    retx[0] = p4x;
    rety[0] = p4y;
    if (count == 2) {
        p4x = p3x - dy;
        p4y = p3y + dx;
        retx[1] = p4x;
        rety[1] = p4y;
    }
    return count;
}

#if (UNITTEST == 1)
int test_circle_get_intersection(){
    int num;
    float ix[2], iy[2];
    
    printf("\nTesting circle_get_intersection\n");
    
    num = circle_get_intersection (0, 0, 1, 0, 1, 1, ix, iy);
    printf("Intersection of (0,0 r:1) (1,0 r:1) should be with %i intersections.",2);
    if (num == 0) printf("found no intersections\n");
    if (num > 0) printf(" - found (%.8f,%.8f)", ix[0], iy[0]);
    if (num == 2) printf("and (%.8f,%.8f)" , ix[1], iy[1]);
    if (num>0) printf(" with %i intersections\n",num);
    if (num!=2) return 0;
    
    num = circle_get_intersection (0, 0, 3, 0, 1, 1, ix, iy);
    printf("Intersection of (0,0 r:1) (3,0 r:1) should be with %i intersections.",0);
    if (num == 0) printf("found no intersections\n");
    if (num > 0) printf(" - found (%.8f,%.8f)", ix[0], iy[0]);
    if (num == 2) printf("and (%.8f,%.8f)" , ix[1], iy[1]);
    if (num>0) printf(" with %i intersections\n",num);
    if (num!=0) return 0;

    num = circle_get_intersection (0, 0, 2, 0, 1, 1, ix, iy);
    printf("Intersection of (0,0 r:1) (1,0 r:1) should be with %i intersections.",1);
    if (num == 0) printf(" - found no intersections\n");
    if (num > 0) printf(" - found (%.8f,%.8f)", ix[0], iy[0]);
    if (num == 2) printf("and (%.8f,%.8f)" , ix[1], iy[1]);
    if (num>0) printf(" with %i intersections\n",num);
    if (num!=1) return 0;
    
    num = circle_get_intersection (400, 400, 650, 410, 408.291016, 658.722534, ix, iy);
    printf("Intersection of (400,400 r:410) (650,410 r:658) should be with %i intersections.",2);
    if (num == 0) printf(" - found no intersections\n");
    if (num > 0) printf(" - found (%.8f,%.8f)", ix[0], iy[0]);
    if (num == 2) printf("and (%.8f,%.8f)" , ix[1], iy[1]);
    if (num>0) printf(" with %i intersections\n",num);
    if (num!=2) return 0;
    
    return 1;
}       
#endif



/**
 * Returns an approximated intersection of the two circles.
 * <p>
 * This method might be helpful if the circles are separate or contained
 * within the other. Returns the intersection of the two circles when
 * equally growing both circles till they intersect in one point.
 *
 * @param p1 The center of the first circle.
 * @param r1 The radius of the first circle.
 * @param p2 The center of the second circle.
 * @param r2 The radius of the second circle.
 *
 * @return An approximated intersection of the two circles or
 *         <code>null</code> if p1 equals p2.
 */
static inline size_t
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
circle_get_approx_intersection(const float  p1x, const float p1y,
                               const float  p2x, const float p2y,
                               const float r1, const float r2,
                               float *restrict retx, float *restrict rety)
{
    // calculate distance between center of circles
    float dist = distance_s(p1x,p1y,p2x,p2y);
    // if distance is zero => infinite number of solutions
    if (dist == 0) return 0;

    // calculate intersection of line through center of circles
    // with both circles => four intersection points
    float dr1 = r1/dist;
    float dr2 = r2/dist;
    float dx = (p2x - p1x);
    float dy = (p2y - p1y);
    float dxp1 = dr1 * dx;
    float dyp1 = dr1 * dy;
    float dxp2 = dr2 * dx;
    float dyp2 = dr2 * dy;

    float p11x = p1x + dxp1;
    float p11y = p1y + dyp1;
    float p12x = p1x - dxp1;
    float p12y = p1y - dyp1;
    float p21x = p2x + dxp2;
    float p21y = p2y + dyp2;
    float p22x = p2x - dxp2;
    float p22y = p2y - dyp2;

    // find nearest pair of intersection points belonging
    // to different circles
    dist = distance_s(p11x,p11y,p21x,p21y);
    float n1x = p11x;
    float n1y = p11y;
    float n2x = p21x;
    float n2y = p21y;
    
    float dt = distance_s(p11x,p11y,p22x,p22y);
    if (dt < dist) {
        dist = dt;
        n2x = p22x;
        n2y = p22y;
    }
    dt = distance_s(p12x,p12y,p21x,p21y);
    if (dt < dist) {
        dist = dt;
        n1x = p12x;
        n1y = p12y;
        n2x = p21x;
        n2y = p21y;
    }
    dt = distance_s(p12x,p12y,p22x,p22y);
    if (dt < dist) {
        n1x = p12x;
        n1y = p12y;
        n2x = p22x;
        n2y = p22y;
    }

    // return middle of line between two nearest points as result
    retx[0] = (n1x + n2x) / 2;
    rety[0] = (n1y + n2y) / 2;
    return 1;
}

    /**
     * Returns an approximated intersection of the two circles as found in
     * paper "A Low-Complexity Geometric Bilateration Method for Localization
     * in Wireless Sensor Networks and Its Comparison with Least-Squares
     * Methods" (Figure 3).
     *
     * @param p1 The center of the first circle.
     * @param r1 The radius of the first circle.
     * @param p2 The center of the second circle.
     * @param r2 The radius of the second circle.
     *
     * @return An approximated intersection of the two circles or
     *         <code>null</code> if p1 equals p2.
     */
static inline size_t
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
circle_get_approx_intersection2(const float  p1x, const float p1y,
                               const float  p2x, const float p2y,
                               const float r1, const float r2,
                               float *restrict retx, float *restrict rety)
{
    // calculate distance between center of circles
    float dist = distance_s(p1x,p1y,p2x,p2y);
    // if distance is zero => infinite number of solutions
    if (dist == 0) return 0;

    // take first circle and calculate new radius for this circle
    // as follows, same for second circle but alleviate because of
    // floating point calculations (use some small epsilon offset)
    // and special case of zero distance:
    double rr1 = r1;
    double rr2 = r2;
    if (rr1 == 0) {
        rr1 = 0.001f;
    }
    if (rr2 == 0) {
        rr2 = 0.001f;
    }
    double r1n = fabs(dist - r2) + 0.001001;
    double r2n = fabs(dist - r1) + 0.001001;
    
    // Use CCI procedure to calculate intersection points
    float ret1x[2];
    float ret1y[2];
    float ret2x[2];
    float ret2y[2];
    memset (ret1x,0,sizeof(float)*2);
    memset (ret1y,0,sizeof(float)*2);
    memset (ret2x,0,sizeof(float)*2);
    memset (ret2y,0,sizeof(float)*2);
    circle_get_intersection (p1x, p1y, p2x, p2y, (float)r1n, (float)rr2, ret1x, ret1y);
    circle_get_intersection (p1x, p1y, p2x, p2y, (float)rr1, (float)r2n, ret2x, ret2y);
    *retx = (ret1x[0] + ret2x[0]) / 2.0f;
    *rety = (ret1y[0] + ret2y[0]) / 2.0f;
    return 1;
}



/**
 * Tests if a given point lies in all circles defined by <code>m</code>
 * and <code>r</code>.
 *
 * @param m The center of the circles.
 * @param r The radius of the circles.
 * @param p The point to be tested.
 *
 * @return <code>true</code> if the given point lies in all circles;
 *         <code>false</code> otherwise.
 */
static inline int
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,__artificial__))
circle_point_in_circles (const int count, const float *restrict mx,
                         const float *restrict my, const float *restrict r,
                         const float px, const float py)
{
    for (int i = 0; i < count; i++) {
        if (distance_s(mx[i],my[i],px,py) > (r[i] + 0.01f))
            return 0; 
    }
    return 1; 
}



/**
 * Same as above but using vector processing. This version uses SSE registers only
 * as geo3 needs it that way. 
 * If the caller wants to test less than 4 circles, it has to pass FLT_MAX into the corresponding
 * maxdists[] fields.
 * 
 * @return <code>true</code> if the given lies inside in all circles;
 *         <code>false</code> otherwise.
 */
static inline int
__attribute__((__always_inline__,__gnu_inline__,__pure__,__artificial__))
circle_point_in_circles_opt128 (__m128 mx, __m128 my, __m128 maxdists,
                                float px, float py)
{
    __m128 vpx = _mm_setr_ps(px, px, px, px);
    __m128 vpy = _mm_setr_ps(py, py, py, py);
#ifdef __AVX__
    return _mm_test_all_ones(_mm_castps_si128(_mm_cmpge_ps(maxdists, distance128(mx, my, vpx, vpy))));        
#else
    return _mm_test_all_ones(_mm_castps_si128(_mm_cmpge_ps(maxdists, distance(mx, my, vpx, vpy))));
#endif
}
 
#endif
