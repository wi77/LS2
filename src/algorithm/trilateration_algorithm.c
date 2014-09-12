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

#ifndef INCLUDED_TRILATERATION_ALGORITHM_C
#define INCLUDED_TRILATERATION_ALGORITHM_C 1

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Trilateration
 ***
 *******************************************************************/
 
 /* @algorithm_name: Trilateration */
 
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
trilateration_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
                  size_t num_anchors __attribute__((__unused__)),
                  size_t width __attribute__((__unused__)),
                  size_t height __attribute__((__unused__)),
                  VECTOR *restrict resx, VECTOR *restrict resy)
{
    /* Fixed Values - but gcc will optimize these calculations out of the loop! */
    const VECTOR pres = vx[2] * vx[2] - vx[1] * vx[1] + vy[2] * vy[2] - vy[1] * vy[1];
    const VECTOR pret = vx[0] * vx[0] - vx[1] * vx[1] + vy[0] * vy[0] - vy[1] * vy[1];
    const VECTOR bigsub = (vx[1] - vx[2]) * (vy[0] - vy[1]) - (vy[2] - vy[1])*(vx[1] - vx[0]);

    VECTOR s = (pres + r[1]*r[1] - r[2]*r[2]) * half;
    VECTOR t = (pret + r[1]*r[1] - r[0]*r[0]) * half;
    *resy = (t * (vx[1] - vx[2])- s * (vx[1] - vx[0])) / bigsub;
    *resx = ((*resy * (vy[0] - vy[1])) - t) / (vx[1] - vx[0]);
}

#endif
