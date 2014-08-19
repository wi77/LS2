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


#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <inttypes.h>
#include <stdbool.h>

#include <glib.h>

#include "ls2/library.h"
#include "ls2/ls2.h"
#include "ls2/backend.h"
#include "ls2/util.h"

#ifndef EPSILON
#  define EPSILON 16
#endif


static char *compare[NUM_VARIANTS];
static double similarity = 3.0;
static double dynamic = 200.0;
double ls2_backend_steps = 0.0;

static GOptionEntry cli_options[] = {
    { "gradation", 'G', 0,
      G_OPTION_ARG_DOUBLE,
      &ls2_backend_steps,
      "number of gradation steps, 0 is unlimited", "steps" },
     { "similarity", 'S', 0, G_OPTION_ARG_DOUBLE,
       &similarity, "similarity threshold", NULL },
     { "dynamic", 'D', 0, G_OPTION_ARG_DOUBLE,
       &dynamic, "dynamic range", NULL },
     { "average", 'o', 0, G_OPTION_ARG_STRING, &(compare[AVERAGE_ERROR]),
       "compare average values.", NULL },
     { "maximum", 'M', 0, G_OPTION_ARG_STRING, &(compare[MAXIMUM_ERROR]),
       "compare maximum values.", NULL },
     { "minimum", 'm', 0, G_OPTION_ARG_STRING, &(compare[MINIMUM_ERROR]),
       "compare minimum values.", NULL },
     { "standard-deviation", 's', 0, G_OPTION_ARG_STRING, &(compare[STANDARD_DEVIATION]),
       "compare standard deviations.", NULL },
     { "rmse", 'r', 0, G_OPTION_ARG_STRING, &(compare[ROOT_MEAN_SQUARED_ERROR]),
       "compare root mean squared errors.", NULL },
     { NULL }
};

int
main(int argc, char **argv)
{
     GOptionContext *opt_con;   /* context for parsing command-line options */
     GError *error;
     uint16_t a_height, a_width, b_height, b_width;
     size_t a_no_anchors, b_no_anchors;
     vector2 *a_anchors, *b_anchors;
     float *a_results, *b_results, *results;

     opt_con = g_option_context_new(" - differences between two spatial "
                                    "distributions");
     g_option_context_add_main_entries(opt_con, cli_options, NULL);

     if (!g_option_context_parse(opt_con, &argc, &argv, &error)) {
          g_print("option parsing failed: %s\n", error->message);
          g_option_context_free(opt_con);
          exit(EXIT_FAILURE);
     }

     if (argc < 3) {
	  fprintf(stderr, "error: missing input file name\n");
	  fflush(stderr);
	  exit(EXIT_FAILURE);
     }

     for (ls2_output_variant var = AVERAGE_ERROR; var < NUM_VARIANTS; var++) {
          if (compare[var] == NULL)
              continue;
	  ls2_hdf5_read_locbased(argv[1], var, &a_anchors, &a_no_anchors,
				 &a_results, &a_width, &a_height);

	  ls2_hdf5_read_locbased(argv[2], var, &b_anchors, &b_no_anchors,
				 &b_results, &b_width, &b_height);
	  if (a_width != b_width || a_height != b_height) {
	       fprintf(stderr, "Sizes differ. Cannot continue.\n");
	       exit(EXIT_FAILURE);
	  }

	  /* Check whether both images are comparable. */
	  if (a_no_anchors != b_no_anchors) {
	       fprintf(stderr, "warning: number of anchors do not match\n");
	  } else {
	       for (size_t i = 0; i < a_no_anchors; i++) {
		    if (a_anchors[i].x != b_anchors[i].x ||
			a_anchors[i].y != b_anchors[i].y) {
			 fprintf(stderr, "warning: anchor positions differ\n");
		    }
	       }
	  }

	  /* Compute the difference between both images.
	   * A positive number means that image one has the larger error.
	   * A negative number means that image two has the larger error.
	   */
	  results = g_new(float, (size_t)(a_width * a_height));
	  for (uint16_t y = 0; y < a_height; y++) {
	       for (uint16_t x = 0; x < a_width; x++) {
		    const int i = x + y * a_width;
		    results[i] = a_results[i] - b_results[i];
	       }
	  }

	  float mu, sigma, min, max;
	  ls2_statistics(results, (size_t) a_width * a_height,
			 &mu, &sigma, &min, &max);
	  fprintf(stdout, "Average difference = %f, sdev = %f, min = %f, "
		  "max = %f\n", mu, sigma, min, max);
	  
	  ls2_cairo_write_png_diff(compare[var], a_anchors, a_no_anchors,
				   results, a_width, a_height,
                                   similarity, dynamic);
	  g_free(results);
     }
     exit(EXIT_SUCCESS);
}
