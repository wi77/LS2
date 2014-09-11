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

#ifdef HAVE_CONFIG_H
#  include <ls2/ls2-config.h>
#endif

#include <glib.h>

#include "crlb_so_algorithm.h"

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Cramer Rao Lower Bound by So
 ***
 *******************************************************************/

/* @algorithm_name: CRLB (H.C. So) */

typedef struct crlb_so_arguments {
	double sdev;
} crlb_so_arguments;


crlb_so_arguments ls2_crlb_so_arguments;


void __attribute__((__nonnull__))
ls2_crlb_so_init_arguments(crlb_so_arguments *arguments)
{
	arguments->sdev = 30.0;
}


GOptionEntry crlb_so_parameters[] = {
	{ "crlb-so-sdev", 0, 0, G_OPTION_ARG_DOUBLE,
          &ls2_crlb_so_arguments.sdev,
          "Standard deviation of the error model", NULL },
        { NULL }
};


void __attribute__((__nonnull__))
ls2_add_crlb_so_option_group(GOptionContext *context)
{
     GOptionGroup *group;
     group = g_option_group_new("crlb-so",
                                "Parameters to the CRLB model of So",
                                "Parameters to the CRLB model of So",
                                NULL, NULL);
     g_option_group_add_entries(group, crlb_so_parameters);
     g_option_context_add_group(context, group);
}




/*!
 * Do not use an error model with this one, it computes nonesense in this
 * situation! Use a constant error of 0.
 *
 * This algorithm calculates the lower bound of the mean-squared
 * estimation error.
 */

#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

/* Compute the Cramer Rao lower bound for TOA based algorithms
 * assuming a normal distribution with standard deviation
 * crlb_so_sdev.
 *
 * The formula is given by H.C. So in Chapter 2 of ...
 */
static inline float __attribute__((__always_inline__,__pure__))
crlb_so_run(const vector2 *anchor, const size_t num_anchors,
            const vector2 *location)
{
    int s;
    float crlb;
    const float v = (float)(ls2_crlb_so_arguments.sdev * ls2_crlb_so_arguments.sdev);
    gsl_matrix *A = gsl_matrix_calloc(2, 2);
    gsl_matrix *Ai = gsl_matrix_alloc(2, 2);
    gsl_permutation *p = gsl_permutation_alloc(2);
    /* See Equation (2.157) of So's Chapter. */
    for (size_t i = 0; i < num_anchors; i++) {
        const float dx = location->x - anchor[i].x;
        const float dy = location->y - anchor[i].y;
        const float d2 = dx * dx + dy * dy;
        A->data[0 /* 0, 0 */] += (dx * dx) / (v * d2);
        A->data[1 /* 1, 0 */] += (dx * dy) / (v * d2);
        A->data[2 /* 0, 1 */] += (dy * dx) / (v * d2);
        A->data[3 /* 1, 1 */] += (dy * dy) / (v * d2);
    }
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_linalg_LU_invert(A, p, Ai);
    /* See Equation (2.158) of So's Chapter. */
    crlb = (float) (gsl_matrix_get(Ai, 0, 0) + gsl_matrix_get(Ai, 1, 1));
    gsl_permutation_free(p);
    gsl_matrix_free(Ai);
    gsl_matrix_free(A);
    return crlb;
}
