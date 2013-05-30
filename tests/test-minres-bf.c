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

#include <stdint.h>

#include <immintrin.h>
#include <stdio.h>
#include <string.h>
#include "vector_shooter.h"
#include "util/util_misc.c"
#include "util/util_vector.c"
#include "algorithm/min_res1_bf_algorithm.h"
#include "algorithm/min_res1_bf_algorithm.c"

static void
test1(void)
{
     const int width = 250;
     const int height = 250;
     const VECTOR ax[4] = { VECTOR_CONST_BROADCAST(50),
			    VECTOR_CONST_BROADCAST(50),
			    VECTOR_CONST_BROADCAST(250),
			    VECTOR_CONST_BROADCAST(200), };
     const VECTOR ay[4] = { VECTOR_CONST_BROADCAST(50),
			    VECTOR_CONST_BROADCAST(200),
			    VECTOR_CONST_BROADCAST(50),
			    VECTOR_CONST_BROADCAST(200), };
     const VECTOR r[4] = { { 106.06f, 126.06f, 156.06f, 106.06f, },
                           { 106.06f, 106.06f, 156.06f, 126.06f, },
                           { 106.06f, 126.06f, 156.06f, 106.06f, },
                           { 106.06f, 106.06f, 156.06f, 126.06f, }, };
     const VECTOR tagx = VECTOR_CONST_BROADCAST(125);
     const VECTOR tagy = VECTOR_CONST_BROADCAST(125);

     VECTOR resx, resy;
     char b1[255], b2[255];

     min_res1_bf_run(ax, ay, r, 4, width, height, &resx, &resy);
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
