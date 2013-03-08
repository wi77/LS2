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
 **  This file is made only for including in the LS² project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/* @algorithm_name: Maximum Likelihood Estimator (Gamma) */

/*******************************************************************
 ***
 ***   Maximum Likelihood Estimator (Gamma);
 ***
 *******************************************************************/

#ifndef MLE_GAMMA_ALGORITHM_C_INCLUDED
#define MLE_GAMMA_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include <assert.h>
#include <gsl/gsl_multimin.h>

#include "algorithm/llsq_algorithm.c"

#ifndef MLE_GAMMA_DEFAULT_SHAPE
#define MLE_GAMMA_DEFAULT_SHAPE 3.0
#endif

#ifndef MLE_GAMMA_DEFAULT_RATE
#define MLE_GAMMA_DEFAULT_RATE (MLE_GAMMA_DEFAULT_SHAPE / 50.0)
#endif

#ifndef MLE_GAMMA_DEFAULT_OFFSET
#define MLE_GAMMA_DEFAULT_OFFSET 0.0
#endif

#ifndef MLE_GAMMA_DEFAULT_STEP
#define MLE_GAMMA_DEFAULT_STEP 25.0
#endif

#ifndef MLE_GAMMA_DEFAULT_EPSILON
#define MLE_GAMMA_DEFAULT_EPSILON 1e-5
#endif

#ifndef MLE_GAMMA_DEFAULT_ITERATIONS
#define MLE_GAMMA_DEFAULT_ITERATIONS 100
#endif

static double mle_gamma_shape      = MLE_GAMMA_DEFAULT_SHAPE;
static double mle_gamma_rate       = MLE_GAMMA_DEFAULT_RATE;
static double mle_gamma_offset     = MLE_GAMMA_DEFAULT_OFFSET;
static double mle_gamma_step       = MLE_GAMMA_DEFAULT_STEP;
static double mle_gamma_epsilon    = MLE_GAMMA_DEFAULT_EPSILON;
static    int mle_gamma_iterations = MLE_GAMMA_DEFAULT_ITERATIONS;


struct poptOption mle_gamma_arguments[] = {
        { "mle-gamma-rate", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_rate, 0,
          "rate of the gamma distribution", NULL },
        { "mle-gamma-shape", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_shape, 0,
          "shape of the gamma distribution", NULL },
        { "mle-gamma-offset", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_offset, 0,
          "offset to the gamma distribution", NULL },
        { "mle-gamma-step", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_step, 0,
          "size of the initial line search step vector", NULL },
        { "mle-gamma-epsilon", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_epsilon, 0,
          "maximum size of the simplex for termination", NULL },
        { "mle-gamma-iterations", 0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gamma_iterations, 0,
          "maximum number of iterations before termination", NULL },
        POPT_TABLEEND
};


struct mle_gamma_point2d {
    double x;
    double y;
};

/* This structure holds the parameters to the likelihood function. */
struct mle_gamma_params {
    struct mle_gamma_point2d *anchors;
    double *ranges;
    size_t no_anchors;
    double factor;
};


static double
__attribute__((__nonnull__))
mle_gamma_likelihood_function(const gsl_vector *X, void *restrict params)
{
    double result = 1.0;
    const struct mle_gamma_params *const p = params;
    for (size_t j = 0; j < p->no_anchors; j++) {
        const double d2 = (gsl_vector_get(X, 0) - p->anchors[j].x) *
                          (gsl_vector_get(X, 0) - p->anchors[j].x) +
                          (gsl_vector_get(X, 1) - p->anchors[j].y) *
                          (gsl_vector_get(X, 1) - p->anchors[j].y);
        const double x = p->ranges[j] + mle_gamma_offset - sqrt(d2);
        if (x < 10.0 * DBL_MIN) {
            result = 0.0;
            break;
        }
        const double probability =
            p->factor * pow(x, mle_gamma_shape - 1.0) * exp(-mle_gamma_rate * x);
        result *= probability;
    }
    return -result;
}

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
mle_gamma_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    /* Step 0: Set up the likelihood function. */
    struct mle_gamma_params p;
    struct mle_gamma_point2d anchors[no_anchors];
    double ranges[no_anchors];
    p.no_anchors = no_anchors;
    p.anchors = anchors;
    for (size_t i = 0; i < no_anchors; i++) {
        p.anchors[i].x = vx[i][0];
        p.anchors[i].y = vy[i][0];
    }
    p.ranges = ranges;
    p.factor = pow(mle_gamma_rate, mle_gamma_shape) / tgamma(mle_gamma_shape);

    gsl_multimin_function f;
    f.f = &mle_gamma_likelihood_function;
    f.n = 2u;
    f.params = &p;

    /* Step 1: Calculate an initial estimate. */
    VECTOR sx, sy;
    llsq_run(vx, vy, r, no_anchors, width, height, &sx, &sy);

    /* Step 2: Call the simplex optimiser. */
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_vector *ss, *x;
    x = gsl_vector_alloc(2);
    ss = gsl_vector_alloc(2);
    
    for (int i = 0; i < VECTOR_OPS; i++) {
        /* Step 2a: Initialize the parameters. */
        for (size_t j = 0; j < no_anchors; j++) {
            p.ranges[j] = r[j][i];
        }
        gsl_vector_set_all(ss, 50.0);
        gsl_vector_set(x, 0, sx[i]);
        gsl_vector_set(x, 1, sy[i]);

        int iter = 0;
        int status;
        double size;

        gsl_multimin_fminimizer_set(s, &f, x, ss);

        /* Step 2b: Iterate the minimization algorithm. */
        do {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);
            if (status)
                break;
            size = gsl_multimin_fminimizer_size(s);
            status = gsl_multimin_test_size(size, mle_gamma_epsilon);
        } while (status == GSL_CONTINUE && iter < mle_gamma_iterations);

        /* Step 2c: Store the result. */
        (*resx)[i] = (float) gsl_vector_get(s->x, 0);
        (*resy)[i] = (float) gsl_vector_get(s->x, 1);
    }

    /* Step 3: Clean up. */
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}

#endif
