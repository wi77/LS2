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
#define _GNU_SOURCE
#endif

#include <immintrin.h>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef HAVE_POPT_H
#  include <popt.h>
#endif

#include "vector_shooter.h"
#include "ls2/library.h"
#include "ls2/ls2.h"
#include "ls2/backend.h"
#include "ls2/util.h"

#include "util/util_math.c"
#include "util/util_vector.c"

#if defined(STAND_ALONE)
#  include INCLUDE_EM_H(ERRORMODEL)
#  include INCLUDE_ALG_H(ALGORITHM)
#endif

#ifndef OUTPUT_DEFAULT
#  define OUTPUT_DEFAULT "result.png"
#endif

#ifndef ALGORITHM_DEFAULT
#  define ALGORITHM_DEFAULT "trilateration"
#endif

#ifndef ERROR_MODEL_DEFAULT
#  define ERROR_MODEL_DEFAULT "eq-noise"
#endif

#ifndef ESTIMATOR_DEFAULT
#  define ESTIMATOR_DEFAULT "gdop"
#endif

#if HAVE_POPT_H
/* Defined in library.c generated by generate.py. */
extern struct poptOption algorithm_arguments[];
extern struct poptOption estimator_arguments[];
extern struct poptOption error_model_arguments[];
#endif

#if defined(ESTIMATOR)
static char const *estimator;
#else
static char const *algorithm;
static char const *error_model;
static int inverted;
static int relative;
static float tag_x;
static float tag_y;
static long seed;
static long runs;
#endif
static int arg_width;
static int arg_height;
static int num_threads;

static char const *output_format;             /* Format of the output files. */ 
static char const *output[NUM_VARIANTS];      /* Names of output files.      */
static char const *output_hdf5;               /* Names of raw output files.  */


int main(int argc, const char* argv[])
{
    struct timeval start_tv, end_tv; /* For measuring the wall-clock time   */
    poptContext opt_con;        /* context for parsing command-line options */
    char const* anchor[MAX_ANCHORS*2+1];/* anchor parameters                */
    float *results[NUM_VARIANTS]; /* Array holding the results. */
#if !defined(ESITMATOR)
    uint64_t *result = NULL;  /* Array holding the result of inverted calculation. */
#endif
    uint16_t no_anchor_args = 0; /* number of anchor parameters seen.       */
#if !defined(ESTIMATOR)
    /* Center of mass of estimations (inverted only) */
    float center_x, sdev_x, center_y, sdev_y;
#endif
    vector2 *anchors;
    int rc;

#if !defined(ESTIMATOR)
    /* Command line options specific to the inverted shooter. */
    static struct poptOption inverted_arguments[] = {
        { "x", 'x', POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &tag_x, 0,
          "X coordinate of the tag.", NULL },
        { "y", 'y', POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &tag_y, 0,
          "Y coordinate of the tag.", NULL },
        POPT_TABLEEND
    };
#endif

    /* Command line arguments. */
    static struct poptOption cli_options[] = {
#if !defined(STAND_ALONE)
#  if !defined(ESTIMATOR)
       { "algorithm", 'a', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &algorithm, 0,
          "selects the algorithm (one of: " ALGORITHMS ")", NULL },
        { "error-model", 'e', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &error_model, 0,
          "selects the error model (one of: " ERROR_MODELS ")", NULL },
        { "inverted", 'i', POPT_ARG_NONE,
          &inverted, 0,
          "run the inverted version of the algorithm", NULL },
        { "relative", 0, POPT_ARG_NONE,
          &relative, 0,
          "Use a relative density representation", NULL },
        { "progress", 'P', POPT_ARG_NONE, &ls2_progress, 0,
          "Periodically report progress", NULL },
#  else
        { "estimator", 'e', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &estimator, 0,
          "selects the error model (one of: " ESTIMATORS ")", NULL },
#  endif
#endif
        { "height", 'h', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &arg_height, 0,
          "height of the playing field.", NULL },
        { "width", 'w', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &arg_width, 0,
          "width of the playing field.", NULL },
#if !defined(ESTIMATOR)
        { "output-average", 'o',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[AVERAGE_ERROR]), 0,
          "name of the average error output image file", "file name" },
        { "output-maximum", 'M',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[MAXIMUM_ERROR]), 0,
          "name of the maximum error output image file", "file name" },
        { "output-minimum", 'm',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[MINIMUM_ERROR]), 0,
          "name of the minimum error output image file", "file name" },
        { "output-sdev", 's',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[STANDARD_DEVIATION]), 0,
          "name of the output image file for the standard deviation",
          "file name" },
        { "output-rmse", 'p',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[ROOT_MEAN_SQUARED_ERROR]), 0,
          "name of the root mean squared error output image", "file name" },
#else
        { "output", 'o',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &(output[ROOT_MEAN_SQUARED_ERROR]), 0,
          "name of the output image file", "file name" },
#endif
        { "output-hdf5", 'H',
          POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
          &output_hdf5, 0,
          "name of the hdf output file for raw result data", "file name" },
#if !defined(ESTIMATOR)
        { "seed", 0, POPT_ARG_LONG, &seed, 0,
          "seed to use for the pseudo random number generators. Default"
          "is based on the current time.",
          "seed" },
        { "runs", 'r', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT,
          &runs, 0,
          "number of runs per pixel (must be divisible by 8)",
          "number of runs" },
#endif
        { "threads", 't', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &num_threads, 0,
          "number of threads to use", "number of threads" },
        { "verbose", 'v', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &ls2_verbose, 0,
          "be verbose, argument indicates the amount of information", NULL },
#if !defined(ESTIMATOR)
        { NULL, '\0', POPT_ARG_INCLUDE_TABLE, inverted_arguments, 0,
          "Options of the inverted algorithm", NULL },
        { NULL, '\0', POPT_ARG_INCLUDE_TABLE, algorithm_arguments, 0,
          "Options of the algorithm", NULL },
        { NULL, '\0', POPT_ARG_INCLUDE_TABLE, error_model_arguments, 0,
          "Options of the error model", NULL },
#else
        { NULL, '\0', POPT_ARG_INCLUDE_TABLE, estimator_arguments, 0,
          "Options of the estimators", NULL },
#endif
        POPT_AUTOHELP
        POPT_TABLEEND
    };

    ls2_verbose = 0;
    arg_width = SIZE;
    arg_height = SIZE;
    num_threads = NUM_THREADS;
    output_format = "png";
#if !defined(ESTIMATOR)
    algorithm = ALGORITHM_DEFAULT;
    error_model = ERROR_MODEL_DEFAULT;
    tag_x = (float) arg_width / 2.0F;
    tag_y = (float) arg_height / 2.0F;
    output[AVERAGE_ERROR] = OUTPUT_DEFAULT;
    runs = RUNS;
    seed = time(NULL);
#else
    estimator = ESTIMATOR_DEFAULT;
    output[ROOT_MEAN_SQUARED_ERROR] = OUTPUT_DEFAULT;
#endif

    opt_con = poptGetContext(NULL, argc, argv, cli_options, 0);
    poptSetOtherOptionHelp(opt_con, "[OPTIONS] [ANCHORX ANCHORY]{3,}");

    // Check for sufficient number of command line arguments
    if (argc < 2) {
        poptPrintUsage(opt_con, stderr, 0);
        poptFreeContext(opt_con);
        exit(EXIT_FAILURE);
    }

    // Parse the command line arguments.
    while ((rc = poptGetNextOpt(opt_con)) >= 0) {
        switch (rc) {
        default:
            break;
        }
    }

    if (rc < -1) {
        /* an error occurred during option processing */
        fprintf(stderr, "%s: %s\n",
                poptBadOption(opt_con, POPT_BADOPTION_NOALIAS),
                poptStrerror(rc));
        poptFreeContext(opt_con);
        exit(EXIT_FAILURE);
    }

    // Handle the left-over arguments.
    while (no_anchor_args < 2 * MAX_ANCHORS && poptPeekArg(opt_con) != NULL) {
        anchor[no_anchor_args] = poptGetArg(opt_con);
        no_anchor_args++;
    }

    if (no_anchor_args < 6) {
        fprintf(stderr, "insufficient number of anchor coordinates.\n");
        poptFreeContext(opt_con);
        exit(EXIT_FAILURE);
    }

    if (no_anchor_args % 2 != 0) {
        fprintf(stderr, "missing anchor coordinate.\n");
        poptFreeContext(opt_con);
        exit(EXIT_FAILURE);
    }

    // parse and normalize anchor coordinates
    const size_t no_anchors = no_anchor_args / 2;

    anchors = calloc(no_anchors, sizeof(vector2));
    for (int i = 0, j = 0; i < no_anchor_args; i += 2, j++) {
	    anchors[j].x = (float) strtoul(anchor[i],   NULL, 10);
	    anchors[j].y = (float) strtoul(anchor[i+1], NULL, 10);
    }

#if !defined(ESTIMATOR)
    int alg = get_algorithm_by_name(algorithm);
    if (alg < 0) {
        fprintf(stderr, "Algorithm \"%s\" unknown, choose one of " ALGORITHMS
                "\n", algorithm);
        exit(EXIT_FAILURE);
    }
    int em = get_error_model_by_name(error_model);
    if (em < 0) {
        fprintf(stderr, "Error model \"%s\" unknown, choose one of "
                ERROR_MODELS "\n", error_model);
        exit(EXIT_FAILURE);
    }
#else
    int est = get_estimator_by_name(estimator);
    if (est < 0) {
        fprintf(stderr, "Estimator \"%s\" unknown, choose one of "
                ESTIMATORS "\n", estimator);
        exit(EXIT_FAILURE);
    }    
#endif

    if (ls2_verbose >= 1) {
#if !defined(ESTIMATOR)
        fprintf(stdout,
                "\nAlgorithm %s, error model %s using %d threads and %ld runs.\n",
		algorithm, error_model, num_threads, runs);
#else
        fprintf(stdout, "\nEstimator %s using %d threads.\n",
		estimator, num_threads);
#endif
        fprintf(stdout, "%zu anchors: (%f; %f)",
                no_anchors, anchors[0].x, anchors[0].y);
        for (size_t i = 1; i < no_anchors; i++) {
            fprintf(stdout, ", (%f; %f)", anchors[i].x, anchors[i].y);
        }
        fprintf(stdout, "\n");
    }

    /* We need at least one thread. */
    if (num_threads < 1)
        num_threads = 1;

    /* Allocate arrays of float forstatistical evaluation. Allocate it
     * if the user requested an image for it or if he wants the data in
     * an HDF5 file. The later case contains all information.
     */
    uint16_t width = (uint16_t) arg_width;
    uint16_t height = (uint16_t) arg_height;
    const size_t sz = ((size_t) width) * ((size_t) height) * sizeof(float);
    memset(results, 0, sizeof(results));
#if !defined(ESTIMATOR)
    if (inverted == 0) {
        for (ls2_output_variant var = 0; var < NUM_VARIANTS; var++) {
            if ((output[var] != NULL && *output[var] != '\0') ||
                (output_hdf5 != NULL && *output_hdf5 != '\0')) {
                if (posix_memalign((void**)&(results[var]), ALIGNMENT, sz) != 0) {
	            perror("posix_memalign()");
	            exit(EXIT_FAILURE);
                }
            }
        }
    } else {
        const size_t s = (size_t) width * (size_t) height * sizeof(uint64_t);
        if (posix_memalign((void**)&(result), ALIGNMENT, s) != 0) {
	    perror("posix_memalign()");
	    exit(EXIT_FAILURE);
        }
    }
#else
    if (posix_memalign((void**)&(results[ROOT_MEAN_SQUARED_ERROR]), ALIGNMENT, sz) != 0) {
      perror("posix_memalign()");
      exit(EXIT_FAILURE);
    }
#endif
    gettimeofday(&start_tv, NULL);

#if !defined(ESTIMATOR)
    /* Sanitize the number of runs. */
    do {
        long t;
        if (inverted == 0)
            t = iceil((long) runs, (long) VECTOR_OPS);
        else
            t = iceil((long) runs, (long) num_threads * VECTOR_OPS);
        if (t != runs) {
    	    runs = t;
    	    fprintf(stderr, "warning: number of runs rounded to %ld\n", runs);
        }
    } while (0);

    if (inverted == 0) {
	if (ls2_progress != 0) {
	    ls2_initialize_progress_bar((size_t) (runs * height * width),
                                        algorithm);
	}
	ls2_distribute_work_shooter(alg, em, num_threads, runs, seed, anchors,
				    no_anchors, results, width, height);
    } else {
	if (ls2_progress != 0) {
	    ls2_initialize_progress_bar((size_t) runs, algorithm);
	}
	ls2_distribute_work_inverted(alg, em, num_threads, runs, seed,
                                     tag_x, tag_y,
				     anchors, no_anchors, result, width,
				     height, &center_x, &sdev_x,
                                     &center_y, &sdev_y);
    }
#else
    ls2_distribute_work_estimator(est, num_threads, anchors, no_anchors,
				  results, width, height);
#endif

    gettimeofday(&end_tv, NULL);

#if !defined(ESTIMATOR)
    if (ls2_progress != 0) {
        ls2_stop_progress_bar();
    }

    // calculate average
    if (inverted == 0) {
	float mu, sigma, min, max;
        if (results[AVERAGE_ERROR] != NULL) {
	    ls2_statistics(results[AVERAGE_ERROR], (size_t) width * height,
		           &mu, &sigma, &min, &max);
	    fprintf(stdout, "MAE = %f, sdev = %f, min = %f, max = %f\n",
		    mu, sigma, min, max);
            fflush(stderr);
        }
    } else {
	fprintf(stdout, "Centroid of location estimations: (%f, %f)"
                        "\n    standard deviations: (%f, %f)\n"
                        "    mean average error: %f\n",
		center_x, center_y, sdev_x, sdev_y,
                distance_s(tag_x, tag_y, center_x, center_y));
    }
#endif

    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    fprintf(stdout, "real %.6f s, user %lu.%06lu s, sys %lu.%06lu s\n",
            (double) (end_tv.tv_sec - start_tv.tv_sec) +
              (double) (end_tv.tv_usec - start_tv.tv_usec) / 1000000.0,
            resources.ru_utime.tv_sec, resources.ru_utime.tv_usec,
            resources.ru_stime.tv_sec, resources.ru_stime.tv_usec);
    fflush(stdout);

#if !defined(ESTIMATOR)
    if (inverted == 0) {
#endif
        for (ls2_output_variant var = 0; var < NUM_VARIANTS; var++) {
            if (output[var] != NULL && *(output[var]) != '\0') {
	      ls2_write_locbased(get_output_format(output_format), output[var],
				 anchors, no_anchors,
				 results[var], width, height);
            }
        }
        if (output_hdf5 != NULL && *output_hdf5 != '\0') {
	    ls2_hdf5_write_locbased(output_hdf5, anchors, no_anchors, results,
                                    width, height);
        }
#if !defined(ESTIMATOR)
    } else {
        if (relative) {
	    ls2_write_inverted(get_output_format(output_format), output[0],
			       0, tag_x, tag_y, anchors, no_anchors,
			       result, width, height, center_x, center_y);
        } else {
	    ls2_write_inverted(get_output_format(output_format), output[0],
			       (uint64_t) runs, tag_x, tag_y, anchors, no_anchors,
			       result, width, height, center_x, center_y);
        }
        if (output_hdf5 != NULL && *output_hdf5 != '\0') {
	    ls2_hdf5_write_inverted(output_hdf5, tag_x, tag_y, anchors,
                                    no_anchors, result, width, height,
                                    center_x, center_y);
        }
    }
#endif
    // clean-ups.
    free(anchors);
    for (ls2_output_variant var = 0; var < NUM_VARIANTS; var++)
        free(results[var]);
    free(result);
    poptFreeContext(opt_con);
    exit(EXIT_SUCCESS);
}
