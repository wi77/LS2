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
 *** Operations on Points. Optimized for SSE.
 ***
 *******************************************************************/

static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
center_of_mass_opt(const __m128 ptsx, const __m128 ptsy, const __m128 mass,
                  float *restrict retx, float *restrict rety)
{
    __m128 m, x, y;
    x = ptsx * mass;
    y = ptsy * mass;
    m = mass;

    x = _mm_hadd_ps(x, x);
    x = _mm_hadd_ps(x, x);
    y = _mm_hadd_ps(y, y);
    y = _mm_hadd_ps(y, y);
    m = _mm_hadd_ps(m, m);
    m = _mm_hadd_ps(m, m);
    
    *retx = x[0]/m[0];
    *rety = y[0]/m[0];
}

static inline int
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
weiszfeld_test_optimum_opt(const __m128 ptsx, const __m128 ptsy,
                           const __m128 weights, const int i)
{   
    __m128 xi = _mm_load1_ps(&ptsx[i]);
    __m128 yi = _mm_load1_ps(&ptsy[i]);
    
#ifdef __AVX__
    __m128 dists = distance128(xi, yi, ptsx, ptsy);
#else
    __m128 dists = distance(xi, yi, ptsx, ptsy);
#endif
    
    __m128 tmpx = weights * ((xi - ptsx) / dists);
    __m128 tmpy = weights * ((yi - ptsy) / dists);
    
    // Remark: division by zero returns infinity, NaN if first parameter is also zero
    // (which is the case for point i itself) or -infinity if first parameter is negative.
    // We make sure these values are set to zero by testing whether they're less than
    // +inf (false for +inf & NaN) and greather than -inf.
    // TODO Flush-To-Zero?
    float tmp = FLT_MAX;
    __m128 fltmax = _mm_load1_ps(&tmp);
    tmpx = _mm_blendv_ps(_mm_setzero_ps(), tmpx, _mm_and_ps(_mm_cmplt_ps(tmpx, fltmax), _mm_cmpgt_ps(tmpx, -fltmax)));
    tmpy = _mm_blendv_ps(_mm_setzero_ps(), tmpy, _mm_and_ps(_mm_cmplt_ps(tmpy, fltmax), _mm_cmpgt_ps(tmpy, -fltmax)));
    
    tmpx = _mm_hadd_ps(tmpx, tmpx);
    tmpx = _mm_hadd_ps(tmpx, tmpx);
    tmpy = _mm_hadd_ps(tmpy, tmpy);
    tmpy = _mm_hadd_ps(tmpy, tmpy);
    
    float result = sqrtf((tmpx[0] * tmpx[0]) + (tmpy[0] * tmpy[0]));
    return (result <= weights[i]) ? i : -1;
}



static inline float
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
weiszfeld_distance_sum_opt(const __m128 mass, const __m128 ptsx,
                           const __m128 ptsy, const float px, const float py,
                           const __m128 epsilon)
{
    
    __m128 vpx = _mm_load1_ps(&px);
    __m128 vpy = _mm_load1_ps(&py);
    __m128 sum = _mm_sqrt_ps((vpx - ptsx) * (vpx - ptsx) + (vpy - ptsy) * (vpy - ptsy) + epsilon);
    sum = _mm_blendv_ps(_mm_setzero_ps(), sum, _mm_cmpgt_ps(mass, _mm_setzero_ps()));
    
    sum = _mm_hadd_ps(sum, sum);
    sum = _mm_hadd_ps(sum, sum);
    
    return sum[0];
}

/**
 * Optimized geometric median function for up two 4 points.
 * If used for 3 points, weights[3] should be zero.
 */
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
point_geometric_median_opt(int count, const __m128 ptsx, const __m128 ptsy,
                           const __m128 weights,
                           float *restrict retx, float *restrict rety)
{
    // Step 1:
    for(int i = 0; i < count; i++) {
        int result = weiszfeld_test_optimum_opt(ptsx,ptsy, weights, i);
        if (result != -1) {
            *retx = ptsx[i];
            *rety = ptsy[i];
            return;
        }
    }
    
    // Step 2:
    float xx,xy;
    center_of_mass_opt(ptsx,ptsy, weights, &xx, &xy);
    
    // Step 3+4:
    float e0, e1;
    float epsilon = 0.000001f;
    __m128 hyperbolaE = _mm_setr_ps(0.001f, 0.001f, 0.001f, 0.001f);
    float xnewx, xnewy;
    int iterations = 0;
    
    e0 = weiszfeld_distance_sum_opt(weights, ptsx,ptsy, xx, xy, hyperbolaE);

#ifdef __AVX__
    __m128 one_sse = _mm_setr_ps(1.0F, 1.0F, 1.0F, 1.0F);
#else
    __m128 one_sse = one;
#endif

    do {
        __m128 vxx = _mm_load1_ps(&xx);
        __m128 vxy = _mm_load1_ps(&xy);
        
#ifdef __AVX__
        __m128 dists = distance128(vxx, vxy, ptsx, ptsy);
#else
        __m128 dists = distance(vxx, vxy, ptsx, ptsy);
#endif
        
        __m128 xt = weights * (ptsx / dists);
        __m128 yt = weights * (ptsy / dists);
        __m128 id = weights * (one_sse / dists);

        xt = _mm_hadd_ps(xt, xt);
        xt = _mm_hadd_ps(xt, xt);
        yt = _mm_hadd_ps(yt, yt);
        yt = _mm_hadd_ps(yt, yt);
        id = _mm_hadd_ps(id, id);
        id = _mm_hadd_ps(id, id);
        
        xnewx = xt[0] / id[0];
        xnewy = yt[0] / id[0];

        e1 = weiszfeld_distance_sum_opt(weights, ptsx,ptsy, xnewx, xnewy, hyperbolaE);
        if (e1 >= e0) break;
        if (((e0 - e1) / e0) < epsilon) break;

        xx = xnewx;
        xy = xnewy;
        e0 = e1;
        iterations++;

    } while (iterations <= 100);

    *retx = xx;
    *rety = xy;
    return;
}
