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

#ifdef LS_WITH_QR
#include <gsl/gsl_linalg.h>
#endif

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

#ifndef MLE_GAMMA_DEFAULT_EPSILON
#define MLE_GAMMA_DEFAULT_EPSILON 1e-2
#endif

#ifndef MLE_GAMMA_DEFAULT_ITERATIONS
#define MLE_GAMMA_DEFAULT_ITERATIONS 100
#endif

static double mle_gamma_shape      = MLE_GAMMA_DEFAULT_SHAPE;
static double mle_gamma_rate       = MLE_GAMMA_DEFAULT_RATE;
static double mle_gamma_offset     = MLE_GAMMA_DEFAULT_OFFSET;
static double mle_gamma_epsilon    = MLE_GAMMA_DEFAULT_EPSILON;
static    int mle_gamma_iterations = MLE_GAMMA_DEFAULT_ITERATIONS;

//#define DEBUG

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
    double gammaval;
    double factor;
};


static double
__attribute__((__nonnull__))
mle_gamma_likelihood_function(const gsl_vector *restrict X, void *restrict params)
{
    double result = 0.0;
    const struct mle_gamma_params *const p = params;
    const double thetaX = gsl_vector_get(X, 0);
    const double thetaY = gsl_vector_get(X, 1);
    for (size_t j = 0; j < p->no_anchors; j++) {
        const double d2 = (thetaX - p->anchors[j].x) *
                          (thetaX - p->anchors[j].x) +
                          (thetaY - p->anchors[j].y) *
                          (thetaY - p->anchors[j].y);
        const double d = sqrt(d2);
        const double Z = p->ranges[j] + mle_gamma_offset - d;
        if (Z <= 0.0) {
            result = INFINITY; // If a measurement is too short, return this.
            break;
        }
        const double likelihood =
            log(p->factor * pow(Z, mle_gamma_shape - 1.0)) - mle_gamma_rate * Z;
        result -= likelihood;
    }
#ifdef DEBUG
    fprintf(stderr, "f(%f, %f) = %f\n", thetaX, thetaY, result);
#endif
    return result;
}



static void
__attribute__((__nonnull__))
mle_gamma_likelihood_gradient(const gsl_vector *restrict X, void *restrict params,
                              gsl_vector *restrict g)
{
    const struct mle_gamma_params *const p = params;
    const double thetaX = gsl_vector_get(X, 0);
    const double thetaY = gsl_vector_get(X, 1);
    double gradX = 0.0, gradY = 0.0;
    for (size_t j = 0; j < p->no_anchors; j++) {
        if (thetaX != p->anchors[j].x && thetaY != p->anchors[j].y) {
            const double d2 = (thetaX - p->anchors[j].x) *
                              (thetaX - p->anchors[j].x) +
                              (thetaY - p->anchors[j].y) *
                              (thetaY - p->anchors[j].y);
            const double d = sqrt(d2);
            const double Z = p->ranges[j] + mle_gamma_offset - d;

            if (__builtin_expect(Z <= 0, 0)) {
                gradX = NAN;
                gradY = NAN;
                break;
            }
            const double t =
                (mle_gamma_rate * (d - Z - mle_gamma_offset) + mle_gamma_shape - 1.0) /
                (d * (d - Z - mle_gamma_offset));
            gradX -= (thetaX -  p->anchors[j].x) * t;
            gradY -= (thetaY -  p->anchors[j].y) * t;
        } // TODO: Otherwise the value of the gradient component is 0?
    }
#ifdef DEBUG
    fprintf(stderr, "df(%f, %f) = (%f, %f)\n", thetaX, thetaY, gradX, gradY);
#endif
    gsl_vector_set(g, 0, gradX);
    gsl_vector_set(g, 1, gradY);
}




static void
__attribute__((__nonnull__,__flatten__,__hot__,__used__))
mle_gamma_likelihood_fdf(const gsl_vector *X, void *restrict params,
                         double *restrict y, gsl_vector *restrict g)
{
    *y = mle_gamma_likelihood_function(X, params);
    mle_gamma_likelihood_gradient(X, params, g);
}



static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
mle_gamma_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors,
              int width __attribute__((__unused__)),
              int height __attribute__((__unused__)),
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
    p.gammaval = tgamma(mle_gamma_shape);
    p.factor = pow(mle_gamma_rate, mle_gamma_shape) / p.gammaval;

    gsl_multimin_function_fdf fdf;
    fdf.n = 2u;
    fdf.f = &mle_gamma_likelihood_function;
    fdf.df = &mle_gamma_likelihood_gradient;
    fdf.fdf = &mle_gamma_likelihood_fdf;
    fdf.params = &p;

    

    /* Step: Call the optimiser. */
    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, 2);

    gsl_vector *x;
    x = gsl_vector_alloc(2);

#if defined(LS2_WITH_QR)
    gsl_vector *b, *tau, *residual;
    b = gsl_vector_alloc(no_anchors - 1);
    tau = gsl_vector_alloc(2);
    residual = gsl_vector_alloc(no_anchors - 1);

    gsl_matrix *A;
    A = gsl_matrix_alloc(no_anchors - 1, 2);
#endif

#if !defined(LS2_WITH_QR)   
    VECTOR sx, sy;
    llsq_run(vx, vy, r, no_anchors, width, height, &sx, &sy);
#endif

    for (int i = 0; i < VECTOR_OPS; i++) {
        /* Step 2a: Initialize the parameters. */
        for (size_t j = 0; j < no_anchors; j++) {
            p.ranges[j] = r[j][i];
        }

#if defined(LS2_WITH_QR)
        /* Step 2b: Calculate an initial estimate.
           We try the linear least squares solution to the inequality
           system defining the domain. */

        for (size_t j = 0; j < no_anchors - 1; j++) {
            gsl_matrix_set(A, j, 0, 2 * (vx[no_anchors - 1][0] - vx[j][0]));
            gsl_matrix_set(A, j, 1, 2 * (vy[no_anchors - 1][0] - vy[j][0]));
            gsl_vector_set(b, j, (r[j][0] + mle_gamma_offset) * (r[j][0] + mle_gamma_offset) - 
                           (r[no_anchors - 1][0] + mle_gamma_offset) * (r[no_anchors - 1][0] + mle_gamma_offset) +
                           vx[no_anchors - 1][0] * vx[no_anchors - 1][0] -
                           vx[j][0] - vx[j][0] +
                           vy[no_anchors - 1][0] * vy[no_anchors - 1][0] -
                           vy[j][0] - vy[j][0]);
        }

        gsl_linalg_QR_decomp(A, tau);
        gsl_linalg_QR_lssolve(A, tau, b, x, residual);
#ifdef DEBUG
        do {
            fprintf(stderr, "x = ");
            gsl_vector_fprintf(stderr, x, "%g");
            fprintf(stderr, "residual = ");
            gsl_vector_fprintf(stderr, residual, "%g");
        } while (0);
#endif
#else
        gsl_vector_set(x, 0, sx[i]);
        gsl_vector_set(x, 1, sy[i]);
#endif

        double likelihood = mle_gamma_likelihood_function(x, &p);
#ifdef DEBUG
        fprintf(stderr, "likelihood = %f\n", likelihood);
#endif
        if (isinf(likelihood))
            goto store_result;

        /* Step 2c: Iterate the minimization algorithm. */
        int iter = 0;
        int status;
        gsl_multimin_fdfminimizer_set(s, &fdf, x, 1e-3, 1e-4);

        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);
            if (status)
                break;
            status = gsl_multimin_test_gradient (s->gradient, mle_gamma_epsilon);
        } while (status == GSL_CONTINUE && iter < mle_gamma_iterations);

        /* Step 2c: Store the result. */
    store_result:
        (*resx)[i] = (float) gsl_vector_get(s->x, 0);
        (*resy)[i] = (float) gsl_vector_get(s->x, 1);
    }

    /* Step 3: Clean up. */
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);

#if defined(LS2_WITH_QR)
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(tau);
    gsl_vector_free(residual);
#endif

}

#endif
