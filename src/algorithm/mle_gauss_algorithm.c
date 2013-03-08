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

/* @algorithm_name: Maximum Likelihood Estimator (Gauss) */

/*******************************************************************
 ***
 ***   Maximum Likelihood Estimator (Gauss);
 ***
 *******************************************************************/

#ifndef MLE_GAUSS_ALGORITHM_C_INCLUDED
#define MLE_GAUSS_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include <assert.h>
#include <gsl/gsl_multimin.h>

#include "algorithm/llsq_algorithm.c"

#ifndef MLE_GAUSS_DEFAULT_MEAN
#define MLE_GAUSS_DEFAULT_MEAN 50.0
#endif

#ifndef MLE_GAUSS_DEFAULT_DEVIATION
#define MLE_GAUSS_DEFAULT_DEVIATION 20.0
#endif

#ifndef MLE_GAUSS_DEFAULT_STEP
#define MLE_GAUSS_DEFAULT_STEP 50.0
#endif

#ifndef MLE_GAUSS_DEFAULT_EPSILON
#define MLE_GAUSS_DEFAULT_EPSILON 1e-2
#endif

#ifndef MLE_GAUSS_DEFAULT_ITERATIONS
#define MLE_GAUSS_DEFAULT_ITERATIONS 200
#endif

static double mle_gauss_mean       = MLE_GAUSS_DEFAULT_MEAN;
static double mle_gauss_deviation  = MLE_GAUSS_DEFAULT_DEVIATION;
static double mle_gauss_step       = MLE_GAUSS_DEFAULT_STEP;
static double mle_gauss_epsilon    = MLE_GAUSS_DEFAULT_EPSILON;
static    int mle_gauss_iterations = MLE_GAUSS_DEFAULT_ITERATIONS;


struct poptOption mle_gauss_arguments[] = {
        { "mle-gauss-deviation", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_deviation, 0,
          "deviation of the gauss distribution", NULL },
        { "mle-gauss-mean", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_mean, 0,
          "mean of the gauss distribution", NULL },
        { "mle-gauss-step", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_step, 0,
          "size of the initial line search step vector", NULL },
        { "mle-gauss-epsilon", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_epsilon, 0,
          "maximum size of the simplex for termination", NULL },
        { "mle-gauss-iterations", 0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_iterations, 0,
          "maximum number of iterations before termination", NULL },
        POPT_TABLEEND
};


struct mle_gauss_point2d {
    double x;
    double y;
};

/* This structure holds the parameters to the likelihood function. */
struct mle_gauss_params {
    struct mle_gauss_point2d *anchors;
    double *ranges;
    size_t no_anchors;
};


static double
__attribute__((__nonnull__))
mle_gauss_likelihood_function(const gsl_vector *X, void *restrict params)
{
    double result = 0.0;
    const struct mle_gauss_params *const p = params;
    const double mean2 = mle_gauss_mean * mle_gauss_mean;
    const double variance = mle_gauss_deviation * mle_gauss_deviation;
    for (size_t j = 0; j < p->no_anchors; j++) {
        const double d2 = (gsl_vector_get(X, 0) - p->anchors[j].x) *
                          (gsl_vector_get(X, 0) - p->anchors[j].x) +
                          (gsl_vector_get(X, 1) - p->anchors[j].y) *
                          (gsl_vector_get(X, 1) - p->anchors[j].y);
        const double d = sqrt(d2);
        const double z = p->ranges[j];
        double logProbability = z * z;
        logProbability -= (2 * z - 2 * mle_gauss_mean) * d;
        logProbability -= 2 * mle_gauss_mean * z;
        logProbability += d2;
        logProbability += variance * (mean2 * log(M_PI) + 2 * log(mle_gauss_deviation) + log(2));
        logProbability /= 2 * variance;
        result += logProbability;
    }
    return result;
}

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
mle_gauss_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    /* Step 0: Set up the likelihood function. */
    struct mle_gauss_params p;
    struct mle_gauss_point2d anchors[no_anchors];
    double ranges[no_anchors];
    p.no_anchors = no_anchors;
    p.anchors = anchors;
    for (size_t i = 0; i < no_anchors; i++) {
        p.anchors[i].x = vx[i][0];
        p.anchors[i].y = vy[i][0];
    }
    p.ranges = ranges;

    gsl_multimin_function f;
    f.f = &mle_gauss_likelihood_function;
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
            status = gsl_multimin_test_size(size, mle_gauss_epsilon);
        } while (status == GSL_CONTINUE && iter < mle_gauss_iterations);

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
