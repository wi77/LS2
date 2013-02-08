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
 *** Misc Helpers
 ***
 *******************************************************************/

// returns the maximum value of all numbers in a vector and a given skalar
static inline float 
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
vector_max_ps(VECTOR x, float max) {
    for (int k = 0; k < VECTOR_OPS; k++)
        if (max < x[k]) max = x[k];
    return max;
}


static inline float 
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
vector_min_ps(VECTOR x, float min)
{
    for (int k = 0; k < VECTOR_OPS; k++)
        if (min > x[k]) min = x[k];
    return min;
}


static inline char*
__attribute__((__always_inline__,__gnu_inline__,__artificial__))
vector_to_nstring(char * restrict buffer, const size_t size,
                  const VECTOR vector)
{
     const char format[] = "(%f, %f, %f, %f"
#ifdef __AVX__
	  ", %f, %f, %f, %f"
#endif
	  ")";

     snprintf(buffer, size, format,
	      vector[0], vector[1], vector[2], vector[3]
#ifdef __AVX__
	      , vector[4], vector[5], vector[6], vector[7]
#endif
	  );
     return buffer;
}


// benchmark helper function 
static inline unsigned long long
__attribute__((__always_inline__,__gnu_inline__,__artificial__))
rdtsc(void)
{
     unsigned a, d;
     asm volatile("cpuid");
     asm volatile("rdtsc" : "=a" (a), "=d" (d));
     return (((unsigned long long)a) | (((unsigned long long)d) << 32));
}
