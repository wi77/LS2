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

/* @algorithm_name: pdop estimation */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct pdop_arguments {
    double scale;
} pdop_arguments;

pdop_arguments ls2_pdop_arguments;

void __attribute__((__nonnull__))
ls2_pdop_init_arguments(struct pdop_arguments *arguments)
{
	arguments->scale = 1.0;
}

GOptionEntry pdop_parameters[] = {
	{ "pdop-scale", 0, 0, G_OPTION_ARG_DOUBLE,
          &ls2_pdop_arguments.scale,
          "scale factor to the root error", NULL },
        { NULL }
};

void __attribute__((__nonnull__))
ls2_add_pdop_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("pdop",
                                "Parameters to the PDOP calculation",
                                "Parameters to the PDOP calculation",
                                NULL, NULL);
     g_option_group_add_entries(group, pdop_parameters);
     g_option_context_add_group(context, group);
}


static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,__artificial__))
pdop_run(const vector2 *restrict anchor, const size_t num_anchors,
         const vector2 *restrict location)
{
    gsl_matrix *A = gsl_matrix_alloc(3, num_anchors);
    for (size_t i = 0; i < num_anchors; i++) {
         const float d = distance_v(location, anchor + i);
         gsl_matrix_set(A, 0, i, (anchor[i].x - location->x) / d);
         gsl_matrix_set(A, 1, i, (anchor[i].y - location->y) / d);
         gsl_matrix_set(A, 2, i, -1.0);
    }

    // transpose
    gsl_matrix *At = gsl_matrix_alloc(num_anchors, 3);
    gsl_matrix_transpose_memcpy(At, A);

    // Multiply
    gsl_matrix *AAt = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, At, 0.0, AAt);

    // invert
    int s;
    gsl_matrix *Q = gsl_matrix_alloc(3, 3);
    gsl_permutation *perm = gsl_permutation_alloc(3);

    gsl_linalg_LU_decomp(AAt, perm, &s);

    // Test if AAt is singular. If so, return NAN
    float pdop;
    for (size_t i = 0; i < 3; i++) {
        if (gsl_matrix_get(AAt, i, i) == 0.0) {
            pdop = NAN;
            goto cleanup;
        }
    }

    gsl_linalg_LU_invert(AAt, perm, Q);
    double trace = 0;
    for (size_t i = 0; i < 2; i++) {
        trace += gsl_matrix_get(Q, i, i);
    }
    pdop = (float) (ls2_pdop_arguments.scale * sqrt(trace));

  cleanup:
    gsl_matrix_free(Q);
    gsl_permutation_free(perm);
    gsl_matrix_free(AAt);
    gsl_matrix_free(At);
    gsl_matrix_free(A);

    return pdop;
}
