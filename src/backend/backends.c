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

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ls2/library.h"
#include "ls2/ls2.h"
#include "ls2/backend.h"

static const char *format_names[] = {
    "png",
    "pdf",
    "OpenEXR",
    NULL
};



enum ls2_output_format_t __attribute__((__const__))
get_output_format(const char *format)
{
    int i;
    for (i = 0; format_names[i] != NULL; i++) {
	if (strcasecmp(format, format_names[i]) == 0) {
		return i;
	}
    }
    fprintf(stderr, "warning: output format \"%s\" unknown (use ", format);
    for (i = 0; format_names[i+1] != NULL; i++)
	fprintf(stderr, "%s, ", format_names[i]);
    fprintf(stderr, "or %s). Using png.\n", format_names[i+1]);
    fflush(stderr);
    return OUTPUT_PNG;
}



void
ls2_write_locbased(ls2_output_format_t format, const char *filename,
		   const vector2 *anchors, const size_t num_anchors,
		   const float *results, const uint16_t width,
		   const uint16_t height)
{
    switch (format) {
    case OUTPUT_PNG:
	ls2_cairo_write_png_locbased(filename, anchors, num_anchors,
				     results, width, height);
	break;
    case OUTPUT_PDF:
	ls2_cairo_write_pdf_locbased(filename, anchors, num_anchors,
				     results, width, height);
	break;
    case OUTPUT_OPENEXR:
	ls2_openexr_write_locbased(filename, anchors, num_anchors,
				   results, width, height);
	break;
    case NUM_OUTPUT_FORMATS:
	break;
    }
}



void
ls2_write_density(ls2_output_format_t format, const char *filename,
		  const vector2 *anchors, const size_t num_anchors,
		  const float *results, const uint16_t width,
		  const uint16_t height)
{
    switch (format) {
    case OUTPUT_PNG:
	ls2_cairo_write_png_density(filename, anchors, num_anchors,
				    results, width, height);
	break;
    case OUTPUT_PDF:
	ls2_cairo_write_pdf_density(filename, anchors, num_anchors,
				    results, width, height);
	break;
    case OUTPUT_OPENEXR:
	// ls2_openexr_write_density(filename, anchors, num_anchors,
	//			  results, width, height);
        abort();
	break;
    case NUM_OUTPUT_FORMATS:
	break;
    }
}




void
ls2_write_inverted(ls2_output_format_t format, const char *filename,
		   const uint64_t runs, const float tag_x, float tag_y,
		   const vector2 *restrict anchors, const size_t num_anchors,
		   const uint64_t *restrict result, const uint16_t width,
		   const uint16_t height,
    		   const float center_x, float center_y)
{
    uint64_t maxval = 0;
    if (runs <= 0) {
        for (size_t i = 0; i < (size_t)(width * height); i++) {
            maxval = MAX(result[i], maxval);
        }
    } else {
        maxval = runs;
    }

    double *converted;
    converted = malloc((size_t)(width * height) * sizeof(*converted));
    if (converted == NULL) {
        perror(__FUNCTION__);
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < (size_t)width * height; i++) {
        if (maxval > 0) {
            double v = ((double) result[i] / (double) maxval);
            converted[i] = CLAMP(0.0, 1.0, v);
        } else {
            converted[i] = 0.0;
        }
    }
    // Switch on format.
    switch (format) {
    case OUTPUT_PNG:
	ls2_cairo_write_png_inverted(filename, tag_x, tag_y, anchors,
				     num_anchors, converted, width, height,
				     center_x, center_y);
	break;
    case OUTPUT_PDF:
	ls2_cairo_write_pdf_inverted(filename, tag_x, tag_y,
				     anchors, num_anchors, converted,
				     width, height, center_x, center_y);
	break;
    case OUTPUT_OPENEXR:
	ls2_openexr_write_inverted(filename, tag_x, tag_y,
				   anchors, num_anchors, converted,
				   width, height, center_x, center_y);
	break;
    case NUM_OUTPUT_FORMATS:
	break;
    }
    free(converted);
}
