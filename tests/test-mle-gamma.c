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

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include <stdint.h>

#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include "vector_shooter.h"
#include "util/util_misc.c"
#include "util/util_vector.c"

#undef NDEBUG
#undef DEBUG_TRACE_GRADIENT
#define DEBUG_TRACE_ITERATIONS
#undef DEBUG_TRACE_LIKELIHOOD

#include "algorithm/mle_gamma_algorithm.h"
#include "algorithm/mle_gamma_algorithm.c"

static void
test1(void)
{
     const int width = 1000;
     const int height = 1000;
#define NO_ANCHORS 3U
     const size_t no_anchors = NO_ANCHORS;
     const VECTOR ax[NO_ANCHORS] = { VECTOR_CONST_BROADCAST(300),
			             VECTOR_CONST_BROADCAST(300),
			             VECTOR_CONST_BROADCAST(700), };
     const VECTOR ay[NO_ANCHORS] = { VECTOR_CONST_BROADCAST(300),
			             VECTOR_CONST_BROADCAST(700),
			             VECTOR_CONST_BROADCAST(300), };
     const VECTOR r[NO_ANCHORS] = { { 305.84f, 282.85f, 282.85f, 282.85f, },
                                    { 486.05f, 282.85f, 282.85f, 282.85f, },
                                    { 80.851f, 282.85f, 282.85f, 282.85f, }, };
     const VECTOR tagx = VECTOR_CONST_BROADCAST(500);
     const VECTOR tagy = VECTOR_CONST_BROADCAST(500);

     VECTOR resx, resy;
     char b1[255], b2[255];

     mle_gamma_run(ax, ay, r, no_anchors, width, height, &resx, &resy);

     printf("resx = %s\nresy = %s\n", vector_to_nstring(b1, 255, resx),
	    vector_to_nstring(b2, 255, resy));
}

int
main(const int argc __attribute__((__unused__)),
     const char *argv __attribute__((__unused__)))
{
     test1();
     return 0;
}
