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
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Constant value localization (debug only)
 ***
 *******************************************************************/

/* @algorithm_name: Constant value localization (debug only) */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
const_run(const VECTOR *restrict vx __attribute__((__unused__)),
          const VECTOR *restrict vy __attribute__((__unused__)),
          const VECTOR *restrict r __attribute__((__unused__)),
          size_t num_anchors __attribute__((__unused__)),
          size_t width, size_t height, VECTOR *restrict resx, VECTOR *restrict resy)
{
    *resx = VECTOR_BROADCASTF((float) width / 2.0f);
    *resy = VECTOR_BROADCASTF((float) height / 2.0f);
}
