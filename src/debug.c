#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#include <immintrin.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include <glib.h>

#include "ls2/library.h"
#include "ls2/ls2.h"

#include "vector_shooter.h"

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"

#include "library.c"

int main(int argc __attribute__((__unused__)),
         char **argv __attribute__((__unused__)))
{
    VECTOR vx[8] = {
	VECTOR_BROADCASTF(100),
	VECTOR_BROADCASTF(100),
	VECTOR_BROADCASTF(400),
	VECTOR_BROADCASTF(800),
	VECTOR_BROADCASTF(110),
	VECTOR_BROADCASTF(110),
	VECTOR_BROADCASTF(410),
	VECTOR_BROADCASTF(810),
    };
    VECTOR vy[8] = {
	VECTOR_BROADCASTF(100),
	VECTOR_BROADCASTF(300),
	VECTOR_BROADCASTF(600),
	VECTOR_BROADCASTF(800),
	VECTOR_BROADCASTF( 90),
	VECTOR_BROADCASTF(290),
	VECTOR_BROADCASTF(590),
	VECTOR_BROADCASTF(790),
    };
    VECTOR r[8] = {
	VECTOR_BROADCASTF(607),
	VECTOR_BROADCASTF(452),
	VECTOR_BROADCASTF(253),
	VECTOR_BROADCASTF(461),
	VECTOR_BROADCASTF(628),
	VECTOR_BROADCASTF(484),
	VECTOR_BROADCASTF(233),
	VECTOR_BROADCASTF(508),
    };
    VECTOR resx, resy;

    if (argc > 2) {
        fprintf(stderr, "Usage: %s [algorithm].\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    int alg = ALG_TRILATERATION; 
    if (argc == 2)
	alg = get_algorithm_by_name(argv[1]);
    if (alg == -1) {
        fprintf(stderr, "Unknown algorithm %s.\nTry one of "
                ALGORITHMS "\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    srand((unsigned int)time(NULL));
    algorithm(alg, vx, vy, r, 8, 1000, 1000, &resx, &resy);
    printf ("\nResult x=%f y=%f\n",resx[1],resy[1]);
}
