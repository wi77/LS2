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


/*******************************************************************
 ***
 *** Operations on Points
 ***
 *******************************************************************/

#ifndef UTIL_POINTS_C_INCLUDED
#define UTIL_POINTS_C_INCLUDED 1

static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
center_of_mass(const int count, const float *restrict ptsx,
              const float *restrict ptsy, const float *restrict mass,
              float *restrict retx, float *restrict rety)
{
    float m = 0, x = 0, y = 0;
    for (int i = 0; i < count; i++) {
        x += ptsx[i] * mass[i];
        y += ptsy[i] * mass[i];
        m += mass[i];
    }
    *retx = x/m;
    *rety = y/m;
}

static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
center_of_mass_eqw(const int count, const float *restrict ptsx,
              const float *restrict ptsy,
              float *restrict retx, float *restrict rety)
{
    float m = 0, x = 0, y = 0;
    for (int i = 0; i < count; i++) {
        x += ptsx[i];
        y += ptsy[i];
        m += 1.0f;
    }
    *retx = x/m;
    *rety = y/m;
}


#if (UNITTEST == 1)
    int test_center_of_mass() {
        printf("\nTesting center_of_mass\n");
        {
        float ptsx[] = {0,3,7};
        float ptsy[] = {0,5,4};
        float weights[] = {1.0,1.0,1.0};
        float retx,rety;
        center_of_mass(3, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        
        {
        float ptsx[] = {0,3,7,17};
        float ptsy[] = {0,5,4,3};
        float weights[] = {1.0,1.0,1.0,1.0};
        float retx,rety;
        center_of_mass(4, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        
        {
        float ptsx[] = {2.5,3.54,7.54,8};
        float ptsy[] = {1.2,5.89,4.9,3.4};
        float weights[] = {1.0,1.0,1.0,1.0};
        float retx,rety;
        center_of_mass(4, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        return 1;
    
    }
#endif 

static inline float
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__pure__,__artificial__))
weiszfeld_distance_sum(const int count, const float *restrict ptsx,
                       const float *restrict ptsy, const float px,
                       const float py, float epsilon)
{
    float sum = 0;
    for (int i = 0; i < count; i++) {
        sum += sqrtf((px - ptsx[i]) * (px - ptsx[i]) + (py - ptsy[i]) * (py - ptsy[i]) + epsilon);
    }
    return sum;
}

static inline int
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__pure__,__artificial__))
weiszfeld_test_optimum(const int count, const float *restrict ptsx,
                       const float *restrict ptsy, const float *restrict weights,
                       const int i)
{
    float sumx = 0;
    float sumy = 0;
    for (int m = 0; m < count; m++) {
        if (m != i) {
            const float dist = distance_s(ptsx[i], ptsy[i], ptsx[m], ptsy[m]);
            if (dist == 0.0f)
                continue;
            sumx += weights[m] * ((ptsx[i] - ptsx[m]) / dist);
            sumy += weights[m] * ((ptsy[i] - ptsy[m]) / dist);
        }
    }
    const float result = sqrtf((sumx * sumx) + (sumy * sumy));
    return (result <= weights[i]) ? i : -1;
}
 

/**
 * Calculate the geometric median of a discrete set of sample points.
 * <p>
 * The geometric median is the point minimizing the sum of distances to the
 * sample points. It is also known as the Fermat–Weber point or 1-median.
 * <p>
 * This method calculates an approximation to the geometric median using
 * Weiszfeld's algorithm.
 *
 * @param pts The set of sample points.
 * @param weights The weight of each point.
 *
 * @return The geometric median of a discrete set of sample points.
 */
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
point_geometric_median(int count, const float *restrict ptsx,
                       const float *restrict ptsy,
                       const float *restrict weights,
                       float *restrict retx, float *restrict rety)
{
    // Step 1:
    for (int i = 0; i < count; i++) {
        int result = weiszfeld_test_optimum(count,ptsx,ptsy, weights, i);
        if (result != -1) {
            *retx = ptsx[i];
            *rety = ptsy[i];
            return;
        }
    }

    // Step 2:
    float xx,xy;
    center_of_mass(count, ptsx,ptsy, weights, &xx, &xy);

    // Step 3+4:
    float e0, e1;
    float epsilon = 0.000001f;
    float hyperbolaE = 0.001f;
    float xnewx, xnewy;
    int iterations = 0;

    do {
        float xt = 0;
        float yt = 0;
        float id = 0;
	float dist = 0;
        for (int i = 0; i < count; i++) {
            dist = distance_s(xx,xy,ptsx[i],ptsy[i]);
            xt += weights[i] * (ptsx[i] / dist);
            yt += weights[i] * (ptsy[i] / dist);
            id += weights[i] * (1.0f / dist);
        }

        xnewx = xt / id;
        xnewy = yt / id;

//	if (isnan(xnewx)) printf("Fuck\n %f %f %i\n",id,dist,count);

        e0 = weiszfeld_distance_sum(count, ptsx,ptsy, xx, xy, hyperbolaE);
        e1 = weiszfeld_distance_sum(count, ptsx,ptsy, xnewx, xnewy, hyperbolaE);
        if (e1 >= e0) break;
        if (((e0 - e1) / e0) < epsilon) break;

        xx = xnewx;
        xy = xnewy;
        iterations++;

    } while (iterations <= 100);
    *retx = xx;
    *rety = xy;
    return;
}

#if (UNITTEST == 1)
    int test_point_geometric_median() {
        printf("\nTesting geometric_median\n");
        {
        float ptsx[] = {0,3,7};
        float ptsy[] = {0,5,4};
        float weights[] = {1.0,1.0,1.0};
        float retx,rety;
        point_geometric_median(3, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        
        {
        float ptsx[] = {0,3,7,17};
        float ptsy[] = {0,5,4,3};
        float weights[] = {1.0,1.0,1.0,1.0};
        float retx,rety;
        point_geometric_median(4, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        
        {
        float ptsx[] = {2.5,3.54,7.54,8};
        float ptsy[] = {1.2,5.89,4.9,3.4};
        float weights[] = {1.0,1.0,1.0,1.0};
        float retx,rety;
        point_geometric_median(4, ptsx, ptsy, weights, &retx, &rety);
        printf ("Median of test case is: %f/%f\n",retx,rety);
        }
        return 1;
    
    }
#endif 

static inline VECTOR 
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
calculate_residual_error(size_t count, const VECTOR *restrict vx,
                         const VECTOR *restrict vy, const VECTOR *restrict r,
                         VECTOR ex, VECTOR ey)
{
    VECTOR error =  VECTOR_BROADCASTF(0.0f);
    for (size_t i = 0; i < count; i++) {
        VECTOR residual = distance(vx[i],vy[i],ex,ey) - r[i];
        error += residual * residual;
    }
    return error;
}

static inline float 
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
calculate_residual_error_s(size_t count, const float *restrict vx,
                         const float *restrict vy, const float *restrict r,
                         float ex, float ey)
{
    float error =  0.0f;
    for (size_t i = 0; i < count; i++) {
        float residual = distance_s(vx[i],vy[i],ex,ey) - r[i];
        error += residual * residual;
    }
    return error;
}
#endif
