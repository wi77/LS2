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

#define _GNU_SOURCE

#include <immintrin.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_POPT_H
#  include <popt.h>
#endif

#include <inttypes.h>

#include "vector_shooter.h"
#include "ls2/ls2.h"
#include "ls2/library.h"

#include "util/util_median.c"
#include "util/util_random.c"
#include "util/util_vector.c"
#include "util/util_matrix.c"
#include "util/util_circle.c"
#include "util/util_vcircle.c"
#include "util/util_triangle.c"
#include "util/util_points.c"
#include "util/util_misc.c"
#include "util/util_points_opt.c"

#include "library.c"

int
main(int argc, const char* argv[])
{
    srand((unsigned int) time(NULL));
    __m128i seed = _mm_set_epi32(rand(), rand(), rand(), rand());
    VECTOR tagx = VECTOR_BROADCASTF(500.0F), tagy = VECTOR_BROADCASTF(500.0F);
    vector2 anchors[4] = { { 300.0, 300.0 }, { 300, 700 }, { 700, 300 }, { 700, 700 } };
    VECTOR vx[4] = { VECTOR_BROADCASTF(300), VECTOR_BROADCASTF(300), VECTOR_BROADCASTF(700), VECTOR_BROADCASTF(700) };
    VECTOR vy[4] = { VECTOR_BROADCASTF(300), VECTOR_BROADCASTF(700), VECTOR_BROADCASTF(300), VECTOR_BROADCASTF(700) };
    VECTOR distances[4];
    VECTOR result[4];
    memset(distances, 0, sizeof(distances));

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <error-model> <samples>.\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    const int em = get_error_model_by_name(argv[1]);
    if (em == -1) {
        fprintf(stderr, "Unknown error model %s.\nTry one of "
                ERROR_MODELS "\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    const long samples = atol(argv[2]);
    if (samples <= 0) {
        fprintf(stderr, "Invalid number of samples %s.\nMust be positive\n",
                argv[2]);
        exit(EXIT_FAILURE);
    }

    error_model_setup(em, anchors, 4);
    for (int i = 0; i < samples; i += VECTOR_OPS) {
	error_model(em, &seed, distances, vx, vy, 4, tagx, tagy, result);
	for (int ii = 0; ii < VECTOR_OPS; ii++) {
	    printf("%f\n", result[0][ii]);
	}
    }
    exit(EXIT_SUCCESS);
}

