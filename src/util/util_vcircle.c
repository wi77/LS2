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

#ifndef UTIL_VCIRCLE_C_INCLUDED
#define UTIL_VCIRCLE_C_INCLUDED 1

/*
 * Compute the instersection of distances, vector version.
 *
 * rx and ry must be two element arrays. The return value is a vector
 * values in {0.0, 1.0, 2.0} indicating the number of intersections.
 * If the result is 0.0 in a component, rx[0..1] and ry[0..1] are both 0.0
 * in the same component. If the result is 1.0 in a component, then rx[0]
 * and ry[0] hold the intersection coordinates in that component, and rx[1]
 * and ry[1] are 0.0 in that component. Otherwise, the components of rx and
 * ry both hold the result coordinates.
 */
static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
vcircle_intersections(const VECTOR p1x, const VECTOR p1y, const VECTOR p2x,
                      const VECTOR p2y, const VECTOR r1, const VECTOR r2,
                      VECTOR *restrict rx, VECTOR *restrict ry)
{
    // Initialise result to zero
    VECTOR result = VECTOR_ZERO();
    rx[0] = rx[1] = ry[0] = ry[1] = VECTOR_ZERO();

    // Compute distances
    const VECTOR d = distance(p1x, p1y, p2x, p2y);

    VECTOR a = (r1 * r1 - r2 * r2 + d * d) / (d + d);
    VECTOR v = r1 * r1 - a * a;
    VECTOR h = VECTOR_SQRT(VECTOR_ABS(v));
    VECTOR dx = (p2x - p1x) / d;
    VECTOR dy = (p2y - p1y) / d;
    VECTOR p3x = p1x + a * dx;
    VECTOR p3y = p1y + a * dy;
    dx *= h;
    dy *= h;
    VECTOR p4x = p3x + dy;
    VECTOR p4y = p3y - dx;

    // Mask vector to mask out nonsensical values. A component is 0x00000000,
    // if the solution is to be ignored and 0xFFFFFFFF if it is valid.
    VECTOR mask = VECTOR_ZERO();

    // Solution, because the circles are not coincident and are contained
    // within each other.

    mask = VECTOR_NE(d, VECTOR_ZERO());
    mask = VECTOR_AND(mask, VECTOR_GE(r1 + r2, d));
    mask = VECTOR_AND(mask, VECTOR_LE(VECTOR_ABS(r1 - r2), d));

    rx[0] = VECTOR_AND(p4x, mask);
    ry[0] = VECTOR_AND(p4y, mask);
    result += VECTOR_AND(one, mask);

    p4x = p3x - dy;
    p4y = p3y + dx;

    // Compute second solution.
    mask = VECTOR_AND(mask, VECTOR_OR(VECTOR_NE(p4x, p3x), VECTOR_NE(p4y, p3y)));
    rx[1] = VECTOR_AND(p4x, mask);
    ry[1] = VECTOR_AND(p4y, mask);
    result += VECTOR_AND(one, mask);

    return result;
}


#if (UNITTEST == 1)
int test_vcircle_intersections(){
    VECTOR num;
    VECTOR ix[2], iy[2];
    float t = 1.0;
    const VECTOR one = VECTOR_BROADCAST(&t);
    t = 2.0;
    const VECTOR two = VECTOR_BROADCAST(&t);
    t = 3.0;
    const VECTOR three = VECTOR_BROADCAST(&t);
    
    printf("\nTesting vcircle_intersections\n");
    
    num = vcircle_intersections(VECTOR_ZERO(), VECTOR_ZERO(), one, VECTOR_ZERO(), one, one, ix, iy);
    printf("Intersection of (0,0 r:1) (1,0 r:1) should have %e intersections.", 2.0);
    if (num[0] == 0.0) printf("found no intersections\n");
    if (num[0] > 0.0) printf(" - found (%.8f,%.8f)", ix[0][0], iy[0][0]);
    if (num[0] == 2.0) printf("and (%.8f,%.8f)" , ix[1][0], iy[1][0]);
    if (num[0] > 0.0) printf(" with %e intersections\n", num[0]);
    if (num[0] != 2.0) return 0;
    
    num = vcircle_intersections(VECTOR_ZERO(), VECTOR_ZERO(), three, VECTOR_ZERO(), one, one, ix, iy);
    printf("Intersection of (0,0 r:1) (3,0 r:1) should have %e intersections.", 0.0);
    if (num[0] == 0.0) printf("found no intersections\n");
    if (num[0] > 0.0) printf(" - found (%.8f,%.8f)", ix[0][0], iy[0][0]);
    if (num[0] == 2.0) printf("and (%.8f,%.8f)" , ix[1][0], iy[1][0]);
    if (num[0] > 0.0) printf(" with %e intersections\n", num[0]);
    if (num[0] != 0.0) return 0;

    num = vcircle_intersections (VECTOR_ZERO(), VECTOR_ZERO(), two, VECTOR_ZERO(), one, one, ix, iy);
    printf("Intersection of (0,0 r:1) (1,0 r:1) should have %e intersections.", 1.0);
    if (num[0] == 0.0) printf(" - found no intersections\n");
    if (num[0] > 0.0) printf(" - found (%.8f,%.8f)", ix[0][0], iy[0][0]);
    if (num[0] == 2.0) printf("and (%.8f,%.8f)" , ix[1][0], iy[1][0]);
    if (num[0] > 0.0) printf(" with %e intersections\n", num[0]);
    if (num[0] != 1.0) return 0;
    
#if 0
    num = vcircle_intersections(400, 400, 650, 410, 408.291016, 658.722534, ix, iy);
    printf("Intersection of (400,400 r:410) (650,410 r:658) should be with %i intersections.",2);
    if (num == 0) printf(" - found no intersections\n");
    if (num > 0) printf(" - found (%.8f,%.8f)", ix[0], iy[0]);
    if (num == 2) printf("and (%.8f,%.8f)" , ix[1], iy[1]);
    if (num>0) printf(" with %i intersections\n",num);
    if (num!=2) return 0;
#endif

    return 1;
}       
#endif

#endif
