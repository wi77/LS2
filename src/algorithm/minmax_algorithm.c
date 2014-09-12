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
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Min-Max aka Bounding Box algorithm.
 ***
 *******************************************************************/

/* @algorithm_name: Min-Max aka Bounding Box algorithm*/

#ifndef LS2_MINMAX_ALGORITHM_C_INCLUDED
#define LS2_MINMAX_ALGORITHM_C_INCLUDED 1

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
minmax_run (const VECTOR *restrict vx, const VECTOR *restrict vy,
            const VECTOR *restrict r, size_t num_anchors,
            size_t width __attribute__((__unused__)),
            size_t height __attribute__((__unused__)),
            VECTOR *restrict resx, VECTOR *restrict resy)
{
    VECTOR east  = vx[0] + r[0];
    VECTOR west  = vx[0] - r[0];
    VECTOR south = vy[0] + r[0];
    VECTOR north = vy[0] - r[0];
    for (size_t i = 1; i < num_anchors; i++) {
        east  = VECTOR_MIN(east,  vx[i] + r[i]);
        west  = VECTOR_MAX(west,  vx[i] - r[i]);
        south = VECTOR_MIN(south, vy[i] + r[i]);
        north = VECTOR_MAX(north, vy[i] - r[i]);
    }
    *resx = east * half + west * half;
    *resy = north * half + south * half;
}

#endif
