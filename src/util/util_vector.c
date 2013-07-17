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
 *** Vector operations
 ***
 *******************************************************************/

#ifndef INCLUDED_UTIL_VECTOR_C
#define INCLUDED_UTIL_VECTOR_C

// Calculates the pair distances of all the points in the 4/8 vectors
static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
distance(const VECTOR vx, const VECTOR vy, const VECTOR wx, const VECTOR wy)
{
   const VECTOR a = vx - wx;
   const VECTOR b = vy - wy;
   return VECTOR_SQRT(a*a + b*b);
}

#ifdef __AVX__
// Calculates the pair distances of all the points in the 4 vectors
// SSE version needed for SSE geometric median.
static inline __m128
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
distance128(const __m128 vx,const __m128 vy, const __m128 wx, const __m128 wy)
{
   const __m128 a = vx - wx;
   const __m128 b = vy - wy;
   return _mm_sqrt_ps(a*a + b*b);
}
#endif

#ifdef VECTOR_MODE_AVX
// Calculates the pair distances of all the points in the 4 vectors
// -- strict SSE version
static inline __m128
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
distance128(const __m128 vx, const __m128 vy, const __m128 wx, const __m128 wy)
{
   __m128 a = vx - wx;
   __m128 b = vy - wy;
   return _mm_sqrt_ps(a*a + b*b);
}
#endif



/*! Calculates the squared distances of two (skalar) points
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__,__nonnull__))
distance_squared_v(const vector2 *v, const vector2 *w)
{
   const float a = v->x - w->x;
   const float b = v->y - w->y;
   return (a * a + b * b);
}




/*! Calculates the distances of two (skalar) points
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__const__,__flatten__,__nonnull__,__artificial__))
distance_v(const vector2 *v, const vector2 *w)
{
   return sqrtf(distance_squared_v(v, w));
}




/*! Calculates the squared distances of two (skalar) points
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
distance_squared_s(const float vx, const float vy, const float wx, const float wy)
{
   const float a = vx - wx;
   const float b = vy - wy;
   return (a * a + b * b);
}




/*! Calculates the distances of two (skalar) points
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__flatten__,__const__,__artificial__))
distance_s(const float vx, const float vy, const float wx, const float wy)
{
   return sqrtf(distance_squared_s(vx, vy, wx, wy));
}




/*! Calculates the squared distances of two (skalar) points
 */
static inline double
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
distance_squared_sf(const double vx, const double vy, const double wx, const double wy)
{
   const double a = vx - wx;
   const double b = vy - wy;
   return (a * a + b * b);
}




/*! Calculates the distances of two (skalar) points
 */
static inline double
__attribute__((__always_inline__,__gnu_inline__,__flatten__,__const__,__artificial__))
distance_sf(const double vx, const double vy, const double wx, const double wy)
{
   return sqrt(distance_squared_sf(vx, vy, wx, wy));
}

#endif
