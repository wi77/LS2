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

#include <glib.h>

#include "ls2/library.h"
#include "ls2/ls2.h"
#include "ls2/backend.h"
#include "ls2/util.h"
#include "vector_shooter.h"

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

/* Defined in library.c generated by generate.py. */
extern GOptionEntry algorithm_arguments[];
extern GOptionEntry estimator_arguments[];
extern GOptionEntry error_model_arguments[];

#if defined(ESTIMATOR)
static char const *estimator;
#else
static char const *algorithm;
static char const *error_model;
static int inverted;
static int relative;
static double tag_x;
static double tag_y;
static long seed;
static long runs;
#endif
static int arg_width;
static int arg_height;
static int num_threads;
#if !defined(ESTIMATOR)
static int ls2_progress;
#endif

static char const *output_format;             /* Format of the output files. */ 
static char const *output[NUM_VARIANTS];      /* Names of output files.      */
static char const *output_hdf5;               /* Names of raw output files.  */

double ls2_backend_steps;



int
main(int argc, char* argv[])
{
    struct timeval start_tv, end_tv; /* For measuring the wall-clock time   */
    GOptionContext *opt_con;    /* context for parsing command-line options */
    GError *error = NULL;
    char const* anchor[MAX_ANCHORS*2+1];/* anchor parameters                */
    float *results[NUM_VARIANTS]; /* Array holding the results. */
#if !defined(ESITMATOR)
    uint64_t *result = NULL;  /* Array holding the result of inverted calculation. */
#endif
    uint16_t no_anchor_args = 0; /* number of anchor parameters seen.       */
#if !defined(ESTIMATOR)
    /* Center of mass of estimations (inverted only) */
    float centre_x, sdev_x, centre_y, sdev_y;
#endif
    vector2 *anchors;

#if !defined(ESTIMATOR)
    /* Command line options specific to the inverted shooter. */
    static GOptionEntry inverted_arguments[] = {
        { "x", 'x', 0, G_OPTION_ARG_DOUBLE, &tag_x,
          "X coordinate of the tag.", NULL },
        { "y", 'y', 0, G_OPTION_ARG_DOUBLE, &tag_y,
          "Y coordinate of the tag.", NULL },
        { NULL }
    };
#endif

    /* Command line arguments. */
    static GOptionEntry entries[] = {
#if !defined(STAND_ALONE)
#  if !defined(ESTIMATOR)
       { "algorithm", 'a', 0, G_OPTION_ARG_STRING, &algorithm,
          "selects the algorithm (one of: " ALGORITHMS ")", NULL },
        { "error-model", 'e', 0, G_OPTION_ARG_STRING, &error_model,
          "selects the error model (one of: " ERROR_MODELS ")", NULL },
        { "inverted", 'i', 0, G_OPTION_ARG_NONE, &inverted,
          "run the inverted version of the algorithm", NULL },
        { "relative", 0, 0, G_OPTION_ARG_NONE, &relative,
          "Use a relative density representation", NULL },
        { "progress", 'P', 0, G_OPTION_ARG_NONE, &ls2_progress,
          "Periodically report progress", NULL },
#  else
        { "estimator", 'e', 0, G_OPTION_ARG_STRING, &estimator,
          "selects the error model (one of: " ESTIMATORS ")", NULL },
#  endif
#endif
        { "height", 'h', 0, G_OPTION_ARG_INT, &arg_height,
          "height of the playing field.", NULL },
        { "width", 'w', 0, G_OPTION_ARG_INT, &arg_width,
          "width of the playing field.", NULL },
        { "gradation", 'G', 0, G_OPTION_ARG_DOUBLE, &ls2_backend_steps,
         "number of gradation steps, 0 is unlimited", "steps" },
#if !defined(ESTIMATOR)
        { "output-average", 'o', 0, G_OPTION_ARG_STRING,
          &(output[AVERAGE_ERROR]),
          "name of the average error output image file", "file name" },
        { "output-maximum", 'M', 0, G_OPTION_ARG_STRING,
          &(output[MAXIMUM_ERROR]),
          "name of the maximum error output image file", "file name" },
        { "output-minimum", 'm', 0, G_OPTION_ARG_STRING,
          &(output[MINIMUM_ERROR]),
          "name of the minimum error output image file", "file name" },
        { "output-sdev", 's', 0, G_OPTION_ARG_STRING,
          &(output[STANDARD_DEVIATION]),
          "name of the output image file for the standard deviation",
          "file name" },
        { "output-rmse", 'p', 0, G_OPTION_ARG_STRING,
          &(output[ROOT_MEAN_SQUARED_ERROR]),
          "name of the root mean squared error output image", "file name" },
#else
        { "output", 'o', 0, G_OPTION_ARG_STRING,
          &(output[ROOT_MEAN_SQUARED_ERROR]),
          "name of the output image file", "file name" },
#endif
        { "output-hdf5", 'H', 0, G_OPTION_ARG_STRING, &output_hdf5,
          "name of the hdf output file for raw result data", "file name" },
#if !defined(ESTIMATOR)
        { "seed", 0, 0, G_OPTION_ARG_INT64, &seed,
          "seed to use for the pseudo random number generators. Default"
          "is based on the current time.",
          "seed" },
        { "runs", 'r', 0, G_OPTION_ARG_INT64, &runs,
          "number of runs per pixel (must be divisible by 8)",
          "number of runs" },
#endif
        { "threads", 't', 0, G_OPTION_ARG_INT, &num_threads,
          "number of threads to use", "number of threads" },
        { "verbose", 'v', 0, G_OPTION_ARG_INT, &ls2_verbose,
          "be verbose, argument indicates the amount of information", NULL },
	{ NULL }
    };

    register_sigsegv_handler();

    ls2_verbose = 0;
    arg_width = SIZE;
    arg_height = SIZE;
#ifdef _SC_NPROCESSORS_ONLN
    num_threads = (int) MIN(1024L, 4 * sysconf(_SC_NPROCESSORS_ONLN));
#else
    num_threads = NUM_THREADS;
#endif
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

    opt_con = g_option_context_new(" - lateration spacial simulator");
    g_option_context_add_main_entries(opt_con, entries, NULL);

#if !defined(ESTIMATOR)
    GOptionGroup *group;
    group = g_option_group_new("inverted",
                               "Set arguments to inverted algorithm",
                               "Set arguments to inverted algorithm",
                               NULL, NULL);
    g_option_group_add_entries(group, inverted_arguments);
    g_option_context_add_group(opt_con, group);

    ls2_add_algorithm_option_groups(opt_con);
    ls2_add_error_model_option_groups(opt_con);
#else
    ls2_add_estimator_option_groups(opt_con);
#endif

    // poptSetOtherOptionHelp(opt_con, "[OPTIONS] [ANCHORX ANCHORY]{3,}");

    if (!g_option_context_parse( opt_con, &argc, &argv, &error)) {
        g_print("option parsing failed: %s\n", error->message);
        g_option_context_free(opt_con);
        exit(EXIT_FAILURE);
    }

    // Handle the left-over arguments.
    no_anchor_args = 0;
    while (no_anchor_args < 2 * MAX_ANCHORS && no_anchor_args < argc) {
        anchor[no_anchor_args] = argv[no_anchor_args];
        no_anchor_args++;
    }

    if (no_anchor_args < 6) {
        g_printerr("insufficient number of anchor coordinates.\n");
        g_option_context_free(opt_con);
        exit(EXIT_FAILURE);
    }

    if (no_anchor_args % 2 == 0) {
        g_printerr("missing anchor coordinate.\n");
        g_option_context_free(opt_con);
        exit(EXIT_FAILURE);
    }

    if (no_anchor_args == 2 * MAX_ANCHORS && no_anchor_args < argc) {
        fprintf(stderr, "too many anchors.\n");
        g_option_context_free(opt_con);
        exit(EXIT_FAILURE);
    }

    // parse and normalize anchor coordinates
    const size_t no_anchors = no_anchor_args / 2;

    anchors = g_new(vector2, no_anchors);
    for (int i = 0, j = 0; i < no_anchor_args; i += 2, j++) {
	    anchors[j].x = (float) strtoul(anchor[i],   NULL, 10);
	    anchors[j].y = (float) strtoul(anchor[i+1], NULL, 10);
    }

#if !defined(ESTIMATOR)
    int alg = get_algorithm_by_name(algorithm);
    if (alg < 0) {
        g_printerr("Algorithm \"%s\" unknown, choose one of " ALGORITHMS
                "\n", algorithm);
        exit(EXIT_FAILURE);
    }
    int em = get_error_model_by_name(error_model);
    if (em < 0) {
        g_printerr("Error model \"%s\" unknown, choose one of "
                ERROR_MODELS "\n", error_model);
        exit(EXIT_FAILURE);
    }
#else
    int est = get_estimator_by_name(estimator);
    if (est < 0) {
        g_printerr("Estimator \"%s\" unknown, choose one of "
                ESTIMATORS "\n", estimator);
        exit(EXIT_FAILURE);
    }    
#endif

    if (ls2_verbose >= 1) {
#if !defined(ESTIMATOR)
        g_print("Algorithm %s, error model %s using %d threads and %ld runs.\n",
		algorithm, error_model, num_threads, runs);
#else
        g_print("Estimator %s using %d threads.\n",
		estimator, num_threads);
#endif
        g_print("%zu anchors: (%f; %f)",
                no_anchors, anchors[0].x, anchors[0].y);
        for (size_t i = 1; i < no_anchors; i++) {
            g_print(", (%f; %f)", anchors[i].x, anchors[i].y);
        }
        g_print("\n");
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
    memset(results, 0, sizeof(results));
#if !defined(ESTIMATOR)
    if (inverted == 0) {
        for (ls2_output_variant var = 0; var < NUM_VARIANTS; var++) {
            if ((output[var] != NULL && *output[var] != '\0') ||
                (output_hdf5 != NULL && *output_hdf5 != '\0')) {
		results[var] = g_new0(float, width * height);
            }
        }
    } else {
	result = g_new(uint64_t, width * height);
    }
#else
    results[ROOT_MEAN_SQUARED_ERROR] = g_new(float, width * height);
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
    	    g_printerr("warning: number of runs rounded to %ld\n", runs);
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
	    char buffer[32];
	    snprintf(buffer, 31, "inverted %s", algorithm);
	    ls2_initialize_progress_bar((size_t) runs, buffer);
	}
	ls2_distribute_work_inverted(alg, em, num_threads, runs, seed,
                                     (float) tag_x, (float) tag_y,
				     anchors, no_anchors, result, width,
				     height, &centre_x, &sdev_x,
                                     &centre_y, &sdev_y);
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
	    g_print("MAE = %f, sdev = %f, min = %f, max = %f\n",
		    mu, sigma, min, max);
        }
    } else {
	g_print("Centroid of location estimations: (%f, %f)"
                "\n    standard deviations: (%f, %f)\n"
                "    mean average error: %f\n",
		centre_x, centre_y, sdev_x, sdev_y,
                distance_s((float) tag_x, (float) tag_y, centre_x, centre_y));
    }
#endif

    struct rusage resources;
    getrusage(RUSAGE_SELF, &resources);
    g_print("real %.6f s, user %lu.%06lu s, sys %lu.%06lu s\n",
            (double) (end_tv.tv_sec - start_tv.tv_sec) +
            (double) (end_tv.tv_usec - start_tv.tv_usec) / 1000000.0,
            resources.ru_utime.tv_sec, resources.ru_utime.tv_usec,
            resources.ru_stime.tv_sec, resources.ru_stime.tv_usec);

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
			       0, (float) tag_x, (float) tag_y,
                               anchors, no_anchors,
			       result, width, height, centre_x, centre_y);
        } else {
	    ls2_write_inverted(get_output_format(output_format), output[0],
			       (uint64_t) runs, (float) tag_x, (float) tag_y,
                               anchors, no_anchors,
			       result, width, height, centre_x, centre_y);
        }
        if (output_hdf5 != NULL && *output_hdf5 != '\0') {
	    ls2_hdf5_write_inverted(output_hdf5, (float) tag_x, (float) tag_y,
				    anchors, no_anchors, result, width, height,
                                    centre_x, centre_y);
        }
    }
#endif
    // clean-ups.
    g_free(anchors);
    for (ls2_output_variant var = 0; var < NUM_VARIANTS; var++)
        g_free(results[var]);
    g_free(result);
    g_option_context_free(opt_con);
    exit(EXIT_SUCCESS);
}
