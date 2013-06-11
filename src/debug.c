#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "vector_shooter.h"
#include "util/util_math.c"
#include "util/util_vector.c"
#include "util/util_circle.c"
#include "util/util_points.c"
#include "algorithm/icla_algorithm.c"

int main() {
    VECTOR vx[4];
    VECTOR vy[4];
    VECTOR r[4];
    VECTOR resx, resy;
    float f;
    // Anchors
    f = 100;
    vx[0] = VECTOR_BROADCAST(&f);
    f = 100;
    vy[0] = VECTOR_BROADCAST(&f);
    f = 100;
    vx[1] = VECTOR_BROADCAST(&f);
    f = 300;
    vy[1] = VECTOR_BROADCAST(&f);
    f = 400;
    vx[2] = VECTOR_BROADCAST(&f);
    f = 600;
    vy[2] = VECTOR_BROADCAST(&f);
    f = 800;
    vx[3] = VECTOR_BROADCAST(&f);
    f = 800;
    vy[3] = VECTOR_BROADCAST(&f);
    // Ranges
    f= 800;
    r[0] = VECTOR_BROADCAST(&f);
    f= 800;
    r[1] = VECTOR_BROADCAST(&f);
    f= 800;
    r[2] = VECTOR_BROADCAST(&f);
    f= 800;
    r[3] = VECTOR_BROADCAST(&f);
    icla_run(vx, vy, r, 4, 1000, 1000, &resx, &resy);
    printf ("\nResult x=%f y=%f\n",resx[1],resy[1]);
}
