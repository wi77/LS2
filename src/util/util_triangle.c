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

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

#ifndef UTIL_TRIANGLE_C_INCLUDED
#define UTIL_TRIANGLE_C_INCLUDED 1

/*******************************************************************
 ***
 *** Triangle operations
 ***
 *******************************************************************/
   
/**
 * Returns the perimeter of a triangle.
 * 
 * @param a Point A of the triangle.
 * @param b Point B of the triangle.
 * @param c Point C of the triangle.
 *
 * @return The perimeter of the triangle.
 */
static inline float __attribute__((__always_inline__,__gnu_inline__,__const__))
triangle_perimeter_s(const float ax, const float ay, const float bx, const float by, const float cx, const float cy)
{
    return distance_s(ax,ay,bx,by) + distance_s(ax,ay,cx,cy) + distance_s(bx,by,cx,cy);
}

#if (UNITTEST == 1)
    int test_triangle_perimeter() {
    float area;
    
    printf("\nTesting triangle_perimeter_s\n");
    
    area = triangle_perimeter_s(0,0,0,4,3,0);
    printf("Area of (0,0) (1,1) (0,1) should be 12.0 - found %.8f\n", area);
    if (fabs(area - 12.0f) > 0.0001){
        printf("fail!\n");
        return 0;
    }
    return 1;
}       
#endif


/*!
 * Returns the area of a triangle.
 *
 * @param a Point A of the triangle.
 * @param b Point B of the triangle.
 * @param c Point C of the triangle.
 *
 * @return The area of the triangle.
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
triangle_area_s(const float ax, const float ay, const float bx,
                const float by, const float cx, const float cy)
{
    // using Heron's formula
    VECTOR v1x = { bx, ax, ax, 0, }; 
    VECTOR v1y = { by, ay, ay, 0, };
    VECTOR v2x = { cx, cx, bx, 0, };
    VECTOR v2y = { cy, cy, by, 0, };
    VECTOR res = distance(v1x, v1y, v2x, v2y);
    res = VECTOR_HADD(res, res);
    res = VECTOR_HADD(res, res);
    float t = 0.5f * VECTOR_GET(res, 0);
    VECTOR s = VECTOR_BROADCASTF(t);
    VECTOR tmp = s - res;
    t = VECTOR_GET(tmp, 3) * VECTOR_GET(tmp, 0) * VECTOR_GET(tmp, 1) *
        VECTOR_GET(tmp, 2);
    s = VECTOR_BROADCASTF(t);
    VECTOR v = VECTOR_SQRT(s);
    return VECTOR_GET(v, 0);
}

#if (UNITTEST == 1)
    int test_triangle_area() {
    float area;
    
    printf("\nTesting triangle_area_s\n");
    
    area = triangle_area_s(0,0,1,1,0,1);
    printf("Area of (0,0) (1,1) (0,1) should be 0.5 - found %.8f\n", area);
    if (fabs(area - 0.5f) > 0.0001){
        printf("fail!\n");
        return 0;
    }
    
    area = triangle_area_s(1,1,4,5,2,3);
    printf("Area of (1,1) (4,5) (2,3) should be 1.0 - found %.8f\n", area);
    if (fabs(area - 1.0f) > 0.0001){
        printf("fail!\n");
        return 0;
    }
    return 1;
}       
#endif

/**
 * Returns the inner circle radius of VECTOR_OPS triangles.
 *
 * @param a Point A of the triangle.
 * @param b Point B of the triangle.
 * @param c Point C of the triangle.
 *
 * @return The inner circle radius.
 */
static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
triangle_icr(const VECTOR ax, const VECTOR ay, const VECTOR bx,
             const VECTOR by, const VECTOR cx, const VECTOR cy)
{
    const VECTOR a = VECTOR_SQRT((bx - ax) * (bx - ax) + (by - ay) * (by - ay));
    const VECTOR b = VECTOR_SQRT((bx - cx) * (bx - cx) + (by - cy) * (by - cy));
    const VECTOR c = VECTOR_SQRT((ax - cx) * (ax - cx) + (ay - cy) * (ay - cy));
    
    const VECTOR perimeter = a + b + c;
    const VECTOR s = perimeter * VECTOR_BROADCASTF(0.5f);
    const VECTOR area = VECTOR_SQRT(s * (s - a) * (s - b) * (s - c));
    
    return area / perimeter;
}

#endif
