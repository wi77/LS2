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

#ifndef UTIL_RESIDUAL_C_INCLUDED
#define UTIL_RESIDUAL_C_INCLUDED 1

static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,
               __artificial__))
residual_s(const float *restrict vx, const float *restrict vy,
           const float *restrict r, size_t no_anchors,
           float resx, float resy)
{
	float residual = 0.0f;
	for (size_t i = 0; i < no_anchors; ++i) {
		float d = distance_s(vx[i], vy[i], resx, resy) - r[i];
		residual += d * d;
        }
	return residual;
}

static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,
               __artificial__))
residual_ps(const VECTOR *restrict vx, const VECTOR *restrict vy,
            const VECTOR *restrict r, size_t no_anchors,
            VECTOR resx, VECTOR resy)
{
	VECTOR residual = VECTOR_ZERO();
	for (size_t i = 0; i < no_anchors; ++i) {
		VECTOR d = distance(vx[i], vy[i], resx, resy) - r[i];
		residual += d * d;
        }
	return residual;
}

static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,
               __artificial__))
residual_vs(int ii, const VECTOR *restrict vx, const VECTOR *restrict vy,
            const VECTOR *restrict r, size_t no_anchors,
            float resx, float resy)
{
	float residual = 0.0f;
	for (size_t i = 0; i < no_anchors; ++i) {
		float d = distance_s(vx[i][ii], vy[i][ii], resx, resy) - r[i][ii];
		residual += d * d;
        }
	return residual;
}

#endif
