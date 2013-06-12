#include <immintrin.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "vector_shooter.h"
#include "util/util_math.c"
#include "util/util_vector.c"
#include "util/util_matrix.c"
#include "util/util_circle.c"
#include "util/util_points.c"
#include "algorithm/llsq_algorithm.c"
#include "algorithm/lms_algorithm.c"
#include <time.h>     

int main() {
    VECTOR vx[8];
    VECTOR vy[8];
    VECTOR r[8];
    VECTOR resx, resy;
    float f;
     srand (time(NULL));
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
    f = 110;
    vx[4] = VECTOR_BROADCAST(&f);
    f = 90;
    vy[4] = VECTOR_BROADCAST(&f);
    f = 110;
    vx[5] = VECTOR_BROADCAST(&f);
    f = 290;
    vy[5] = VECTOR_BROADCAST(&f);
    f = 410;
    vx[6] = VECTOR_BROADCAST(&f);
    f = 590;
    vy[6] = VECTOR_BROADCAST(&f);
    f = 810;
    vx[7] = VECTOR_BROADCAST(&f);
    f = 790;
    vy[7] = VECTOR_BROADCAST(&f);
    // Ranges
    f= 607;
    r[0] = VECTOR_BROADCAST(&f);
    f= 452;
    r[1] = VECTOR_BROADCAST(&f);
    f= 253;
    r[2] = VECTOR_BROADCAST(&f);
    f= 461;
    r[3] = VECTOR_BROADCAST(&f);
    f= 628;
    r[4] = VECTOR_BROADCAST(&f);
    f= 484;
    r[5] = VECTOR_BROADCAST(&f);
    f= 233;
    r[6] = VECTOR_BROADCAST(&f);
    f= 508;
    r[7] = VECTOR_BROADCAST(&f);
    lms_run(vx, vy, r, 8, 1000, 1000, &resx, &resy);
    printf ("\nResult x=%f y=%f\n",resx[1],resy[1]);
}
