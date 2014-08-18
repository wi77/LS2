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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define _XOPEN_SOURCE 600

#include <assert.h>
#include <stdbool.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <OpenEXR/ImfCRgbaFile.h>

#include "ls2/library.h"
#include "ls2/ls2.h"
#include "vector_shooter.h"
#include "ls2/backend.h"

#include "backend/colors.c"
#include "util/util_colors.c"

static inline void
__attribute__((__gnu_inline__,__always_inline__,__artificial__))
paint_dot_to_openexr(ImfRgba *image, uint16_t width, uint16_t height,
                     int x, int y, float r, float g, float b)
{
    if (0 <= x && x < width && 0 <= y && y < height) {
        ImfRgba *dot = &(image[x + y * width]);
        ImfFloatToHalf(r, &(dot->r));
        ImfFloatToHalf(g, &(dot->g));
        ImfFloatToHalf(b, &(dot->b));
        ImfFloatToHalf(1.0F, &(dot->a));
   }
}


#define DRAW_PIXEL(x, y) \
  paint_dot_to_openexr(image, width, height, x, y, r, g, b);

#include "bits/circle.c"

static
DECLARE_DRAW_CIRCLE(ls2_draw_circle_to_openexr, ImfRgba *, float r, float g, float b)

#include "bits/rectangle.c"

static
DECLARE_DRAW_RECTANGLE(ls2_draw_rectangle_to_openexr, ImfRgba *, float r, float g, float b)


static void __attribute__((__nonnull__))
ls2_draw_result_to_openexr(ImfRgba *restrict image,
                           const uint16_t width, const uint16_t height,
                           const float *restrict result)
{
    for (uint16_t y = 0; y < height; y++) {
        for (uint16_t x = 0; x < width; x++) {
            const int pos = x + width * y;
	    double r, g, b, a;
            const float sample = result[pos];
	    ls2_pick_color_locbased(sample, &r, &g, &b, &a);
	    paint_dot_to_openexr(image, width, height, x, y,
                                 (float) r, (float) g, (float) b);
        }
    }
}


static void
ls2_draw_inverted_to_openexr(ImfRgba *image,
			     const uint16_t width, const uint16_t height,
			     const double* result)
{
    for (uint16_t y = 0; y < height; y++) {
        for (uint16_t x = 0; x < width; x++) {
            const int pos = x + width * y;
	    double h, s, l;
	    float r, g, b;
	    ls2_pick_color_inverted((double) result[pos], &h, &s, &l);
	    hsl_to_rgbf((float)h, (float)s, (float)l, &r, &g, &b);
	    paint_dot_to_openexr(image, width, height, x, y, r, g, b);
        }
    }
}


static void
ls2_draw_anchors_to_openexr(ImfRgba *image,
                            const uint16_t width, const uint16_t height,
                            const vector2 *anchors, const size_t no_anchors, 
                            const uint16_t offset)
{
    for (size_t i = 0; i < no_anchors; i++) {
	ls2_draw_rectangle_to_openexr(image, width, height,
				      (uint16_t) (anchors[i].x + offset),
				      (uint16_t) anchors[i].y, 2, 2,
				      1.0F, 0.0F, 0.0F);
    }
}


extern void
ls2_openexr_write_locbased(const char* filename,
			   const vector2 *anchors, const size_t no_anchors, 
			   const float* restrict result, const uint16_t width,
			   const uint16_t height)
{
    ImfRgba *image = (ImfRgba *) malloc((size_t) width * height * sizeof(ImfRgba));
    if (image == NULL) {
        perror("ls2_write_openexr(): malloc()");
        exit(EXIT_FAILURE);
    }

    ls2_draw_result_to_openexr(image, width, height, result);
    ls2_draw_anchors_to_openexr(image, width, height, anchors, no_anchors, 0);

    ImfHeader *header = ImfNewHeader();
    ImfOutputFile *file = ImfOpenOutputFile(filename, header, IMF_WRITE_RGBA);
    ImfOutputSetFrameBuffer (file, image, 1, width);
    ImfOutputWritePixels (file, height);

  // finish:
    ImfCloseOutputFile(file);
    ImfDeleteHeader(header);
}


extern void __attribute__((nonnull))
ls2_openexr_write_inverted(const char* filename,
                           const float tag_x, const float tag_y,
			   const vector2 *anchors, const size_t no_anchors,
			   const double *restrict result,
			   const uint16_t width, const uint16_t height,
    			   const double center_x, const double center_y)
{
    ImfRgba *image = (ImfRgba *) calloc((size_t) width * height, sizeof(ImfRgba));
    if (image == NULL) {
        perror("ls2_write_inverted_openexr(): calloc()");
        exit(EXIT_FAILURE);
    }

    // Draw the result image.
    ls2_draw_inverted_to_openexr(image, width, height, result);
    ls2_draw_circle_to_openexr(image, width, height, (uint16_t) tag_x, (uint16_t) tag_y,
                               50, 0.0F, 1.0F, 0.0F);
    ls2_draw_rectangle_to_openexr(image, width, height, (uint16_t) tag_x,
                                  (uint16_t) tag_y, 2, 2, 1.0F, 1.0F, 0.0F);
    ls2_draw_anchors_to_openexr(image, width, height, anchors, no_anchors, 0);
    if (0.0 <= center_x && center_x < width &&
        0.0 <= center_y && center_y < height) {
        ls2_draw_rectangle_to_openexr(image, width, height,
                                      (uint16_t) center_x, (uint16_t) center_y,
				      2, 2, 1.0F, 0.0F, 1.0F);
    }

    ImfHeader *header = ImfNewHeader();
    ImfHeaderSetDisplayWindow(header, 0, 0, width - 1, height - 1);
    ImfHeaderSetDataWindow(header, 0, 0, width - 1, height - 1);
    ImfHeaderSetScreenWindowWidth(header, (float) width - 1);

    ImfOutputFile *file = ImfOpenOutputFile(filename, header, IMF_WRITE_RGBA);
    ImfOutputSetFrameBuffer(file, image, 1, width);
    ImfOutputWritePixels(file, height);

  // finish:
    ImfCloseOutputFile(file);
    ImfDeleteHeader(header);
    free(image);
}
