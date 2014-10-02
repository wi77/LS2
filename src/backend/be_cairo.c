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

#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>

#include <glib.h>

#include <cairo.h>
#include <cairo/cairo-pdf.h>

#include "ls2/library.h"
#include "ls2/output.h"
#include "vector_shooter.h"
#include "ls2/backend.h"

#include "backend/colors.c"
#include "util/util_colors.c"



/*!
 * Draw a result image to a surface.
 */
static cairo_surface_t *
ls2_draw_result_to_cairo(cairo_surface_t *surface,
		         const float *result, const uint16_t width,
			 const uint16_t height,
                         const float expectation, const float clip)
{
    cairo_t *cr;

    // Create a drawing buffer.
    cr = cairo_create(surface);
    
    // Fill the background
    cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    // Color each location by the average or maximum error
    for (uint16_t y = 0; y < height; y++) {
	for (uint16_t x = 0; x < width; x++) {
	    const float sample = result[y * width + x];
            double r, g, b, a;
            ls2_pick_color_locbased(&r, &g, &b, &a, sample, expectation, clip);
	    cairo_set_source_rgb(cr, r, g, b);
	    cairo_rectangle(cr, x, y, 1.0, 1.0);
	    cairo_fill(cr);
	}
    }

    cairo_destroy(cr);
    cairo_surface_flush(surface);

    return surface;
}


static void
ls2_draw_anchors_to_cairo(cairo_surface_t *surface,
		          const vector2 *anchors, const size_t num_anchors, 
                          const uint16_t offset, const double font_size)
{
    cairo_t *cr;
    char buffer[256];

    cr = cairo_create(surface);

    cairo_new_path(cr);
    cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
    cairo_select_font_face(cr, "sans-serif",
			   CAIRO_FONT_SLANT_NORMAL,
			   CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, font_size);
    for (size_t i = 0; i < num_anchors; ++i) {
	// Draw a circle at the anchors position.
	const double x = anchors[i].x + offset;
	const double y = anchors[i].y;
	cairo_new_path(cr);
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	cairo_arc(cr, x, y, 3.0, 0, 2.0 * M_PI);
        cairo_fill(cr);
	cairo_stroke(cr);
	
	// Label the circle
	cairo_text_extents_t te;
	cairo_new_path(cr);
	cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
	snprintf(buffer, 8, "%zu", i + 1);

	cairo_text_extents(cr, buffer, &te);
	cairo_move_to(cr, x - te.width - te.x_bearing - 8.0,
		      y - te.height - te.y_bearing - 8.0);
	cairo_show_text(cr, buffer);
    }
    cairo_surface_flush(surface);
    cairo_destroy(cr);
}


static void
ls2_draw_to_cairo_surface_locbased(cairo_surface_t *surface,
				   const vector2 *anchors, const size_t no_anchors,
				   const float* result, const uint16_t width,
				   const uint16_t height,
                                   const float expectation, const float clip)
{
    const double fn_size = (double) ((width < height) ? width : height) / 50.0;
    ls2_draw_result_to_cairo(surface, result, width, height, expectation, clip);
    ls2_draw_anchors_to_cairo(surface, anchors, no_anchors, 0, fn_size);
                              
}


extern void
ls2_cairo_write_pdf_locbased(const char* filename,
		             const vector2 *anchors, const size_t no_anchors, 
		             const float* result, const uint16_t width,
		             const uint16_t height,
                             const float expectation, const float clip)
{
    cairo_surface_t *surface =
        cairo_pdf_surface_create(filename, width, height);
    ls2_draw_to_cairo_surface_locbased(surface, anchors, no_anchors, result,
				       width, height, expectation, clip);
    cairo_surface_destroy(surface);
}


extern void
ls2_cairo_write_png_locbased(const char* filename,
		             const vector2 *anchors, const size_t no_anchors, 
		             const float* result, const uint16_t width,
		             const uint16_t height,
                             const float expectation, const float clip)
{
    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    ls2_draw_to_cairo_surface_locbased(surface, anchors, no_anchors, result,
				       width, height, expectation, clip);
    cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
}



/*!
 * Draw a result image to a surface.
 */
static cairo_surface_t *
ls2_draw_density_to_cairo(cairo_surface_t *surface,
		         const float *result, const uint16_t width,
			 const uint16_t height)
{
    cairo_t *cr;

    // Create a drawing buffer.
    cr = cairo_create(surface);
    
    // Fill the background
    cairo_set_source_rgb(cr, 0.0, 0.0, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    // Color each location
    for (uint16_t y = 0; y < height; y++) {
	for (uint16_t x = 0; x < width; x++) {
	    const float sample = result[y * width + x];
            double hue, lightness, saturation, r, g, b;
            ls2_pick_color_density(sample, &hue, &saturation, &lightness);
            hsl_to_rgb(hue, saturation, lightness, &r, &g, &b);
	    cairo_set_source_rgb(cr, r, g, b);
	    cairo_rectangle(cr, x, y, 1.0, 1.0);
	    cairo_fill(cr);
	}
    }

    cairo_destroy(cr);
    cairo_surface_flush(surface);

    return surface;
}


static void
ls2_draw_to_cairo_surface_density(cairo_surface_t *surface,
				  const vector2 *anchors, const size_t no_anchors,
				  const float* result, const uint16_t width,
				  const uint16_t height)
{
    const double fn_size = (double) ((width < height) ? width : height) / 50.0;
    ls2_draw_density_to_cairo(surface, result, width, height);
    ls2_draw_anchors_to_cairo(surface, anchors, no_anchors, 0, fn_size);
                              
}


extern void
ls2_cairo_write_pdf_density(const char* filename,
			    const vector2 *anchors, const size_t no_anchors, 
			    const float* result, const uint16_t width,
			    const uint16_t height)
{
    cairo_surface_t *surface =
        cairo_pdf_surface_create(filename, width, height);
    ls2_draw_to_cairo_surface_density(surface, anchors, no_anchors, result,
				      width, height);
    cairo_surface_destroy(surface);
}





extern void
ls2_cairo_write_png_density(const char* filename,
			    const vector2 *anchors, const size_t no_anchors, 
			    const float* result, const uint16_t width,
			    const uint16_t height)
{
    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    ls2_draw_to_cairo_surface_density(surface, anchors, no_anchors, result,
				      width, height);
    cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
}





/*
 * Draw an arrow.
 */
static void __attribute__((__nonnull__))
ls2_cairo_draw_arrow(cairo_t *cr, double x, double y, double end_x, double end_y, double r, double g, double b)
{
	cairo_set_source_rgb(cr, r, g, b);
	cairo_arc(cr, x, y, 1.5, 0, 2.0 * M_PI);
        cairo_fill(cr);
	cairo_stroke(cr);

	const double angle = atan2(end_x - x, end_y - y) + M_PI;
	const double head_length = 5.0;
	const double head_angle = M_PI * 30.0 / 180.0;

	cairo_move_to(cr, x, y);
	cairo_line_to(cr, end_x, end_y);
	cairo_stroke(cr);
	cairo_move_to(cr, end_x, end_y);
	cairo_rel_line_to(cr, head_length * cos(angle - head_angle),
			  head_length * sin(angle - head_angle));
	cairo_stroke(cr);
	cairo_move_to(cr, end_x, end_y);
	cairo_rel_line_to(cr, head_length * cos(angle + head_angle),
			  head_length * sin(angle + head_angle));
	cairo_stroke(cr);
}




/*
 *
 */
struct ls2_phase_portrait_item {
    float length;
    double sx, sy, ex, ey;
};



static int
ls2_phase_portrait_item_cmp(const void *__x, const void *__y)
{
    const struct ls2_phase_portrait_item* x = __x;
    const struct ls2_phase_portrait_item* y = __y;
    if (x->length > y->length)
        return -1;
    else if (x->length < y->length)
        return 1;
    else
        return 0;
}




/*!
 * Draw a result image to a surface.
 */
static cairo_surface_t *
ls2_cairo_draw_phase_portrait(cairo_surface_t *surface,
			      const vector2 *anchors,
			      const size_t no_anchors, 
			      const float *dx, const float *dy,
			      const uint16_t width, const uint16_t height,
			      const uint16_t stride, const float expectation,
                              const float clip)
{
    cairo_t *cr;
    size_t p = 0;
    const size_t no_items =
        (size_t) (((width + 1U) / stride) * ((height + 1U) / stride));
    struct ls2_phase_portrait_item items[no_items];

    // Create a drawing buffer.
    cr = cairo_create(surface);
    
    // Fill the background
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    /* Color each location by the average or maximum error
     *
     * Draw largest arrows first, and finish with drawing the
     * shortest last.
     */
    for (int y = stride / 2; y < height; y += stride) {
	for (int x = stride / 2; x < width; x += stride) {
            int pos = x + y * width;

            items[p].length = sqrtf(dx[pos] * dx[pos] + dy[pos] * dy[pos]);
            items[p].sx     = (double) x;
            items[p].sy     = (double) y;
            items[p].ex     = (double) x + dx[pos];
            items[p].ey     = (double) y + dy[pos];
            p++;
        }
    }
    qsort(items, no_items, sizeof(struct ls2_phase_portrait_item),
          ls2_phase_portrait_item_cmp);

    for (size_t i = 0; i < no_items; i++) {
        double r, g, b, a;

        ls2_pick_color_locbased(&r, &g, &b, &a, items[i].length,
                                expectation, clip);
        ls2_cairo_draw_arrow(cr, items[i].sx, items[i].sy,
                             items[i].ex, items[i].ey, r, g, b);
    }
    cairo_destroy(cr);
    cairo_surface_flush(surface);
    const double fn_size = (double) ((width < height) ? width : height) / 50.0;
    ls2_draw_anchors_to_cairo(surface, anchors, no_anchors, 0, fn_size);
    return surface;
}


void
ls2_cairo_write_png_phase_portrait(const char* filename,
				   const vector2 *anchors,
				   const size_t no_anchors, 
				   const float* dx, const float* dy,
				   const uint16_t width, const uint16_t height,
				   const uint16_t stride,
                                   const float expectation, const float clip)
{
    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    ls2_cairo_draw_phase_portrait(surface, anchors, no_anchors, dx, dy,
				  width, height, stride, expectation, clip);
    cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
}



void
ls2_cairo_write_pdf_phase_portrait(const char* filename,
				   const vector2 *anchors,
				   const size_t no_anchors, 
				   const float* dx, const float* dy,
				   const uint16_t width, const uint16_t height,
				   const uint16_t stride,
                                   const float expectation, const float clip)
{
    cairo_surface_t *surface =
        cairo_pdf_surface_create(filename, width, height);
    ls2_cairo_draw_phase_portrait(surface, anchors, no_anchors, dx, dy,
				  width, height, stride, expectation, clip);
    cairo_surface_destroy(surface);
}



static void
ls2_cairo_draw_diff(cairo_surface_t *surface,
		    const vector2 *anchors, const size_t no_anchors,
                    const float* result, const uint16_t width,
		    const uint16_t height, const double similar,
                    const double dynamic)
{
    cairo_t *cr;

    // Create a drawing buffer.
    cr = cairo_create(surface);
    
    // Fill the background
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    for (uint16_t y = 0; y < height; y++) {
	for (uint16_t x = 0; x < width; x++) {
            double r, g, b, lightness, saturation, hue;
	    const float sample = result[y * width + x];
	    ls2_pick_color_diff(sample, similar, dynamic,
				&hue, &saturation, &lightness);
            hsl_to_rgb(hue, saturation, lightness, &r, &g, &b);
	    cairo_set_source_rgb(cr, r, g, b);
	    cairo_rectangle(cr, x, y, 1.0, 1.0);
	    cairo_fill(cr);
	}
    }

    cairo_destroy(cr);
    cairo_surface_flush(surface);

    const double fn_size = (double) ((width < height) ? width : height) / 50.0;
    ls2_draw_anchors_to_cairo(surface, anchors, no_anchors, 0, fn_size);
}


void
ls2_cairo_write_png_diff(const char* filename,
		         const vector2 *anchors, const size_t no_anchors, 
		         const float* result, const uint16_t width,
		         const uint16_t height,
                         const float similar, const float dynamic)
{
    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    ls2_cairo_draw_diff(surface, anchors, no_anchors, result, width, height,
                        similar, dynamic);
    cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
}



/*!
 * Draw a result image of an inverted computation to a surface.
 */
static void __attribute__((__nonnull__,__flatten__))
ls2_draw_inverted_result_to_cairo(cairo_surface_t *surface,
				  const double *restrict result,
				  const uint16_t width, const uint16_t height)
{
    cairo_t *cr;

    // Create a drawing buffer.
    cr = cairo_create(surface);
    
    // Fill the background
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    // Color each location by the average or maximum error
    for (uint16_t y = 0; y < height; y++) {
	for (uint16_t x = 0; x < width; x++) {
            double h, s, l, r, g, b;
	    const double sample = result[y * width + x];
            if (sample > 0.0) {
                ls2_pick_color_inverted(sample, &h, &s, &l);
                hsl_to_rgb(h, s, l, &r, &g, &b);
		cairo_set_source_rgb(cr, r, g, b);
	        cairo_rectangle(cr, (double) x, (double) y, 1.0, 1.0);
	        cairo_fill(cr);
            }
	}
    }

    // We finished writing the image.
    cairo_destroy(cr);
    cairo_surface_flush(surface);
}


static void __attribute__((__nonnull__))
ls2_inverted_to_cairo_surface(cairo_surface_t *restrict surface,
                              const float tag_x, const float tag_y,
			      const vector2 *restrict anchors,
                              const size_t no_anchors, 
			      const double *restrict result,
			      const uint16_t width, const uint16_t height,
			      const double center_x, const double center_y,
                              const double expectation,
                              const double clip __attribute__((__unused__)))
{
    ls2_draw_inverted_result_to_cairo(surface, result, width, height);
    cairo_t *cr = cairo_create(surface);

    // Add a yellow circle into the center
    cairo_new_path(cr);
    cairo_set_line_width(cr, 2.0);
    cairo_set_source_rgb(cr, 1.0, 1.0, 0.0);
    cairo_arc(cr, tag_x, tag_y, 3.0, 0.0, 2.0 * M_PI);
    cairo_stroke(cr);
    
    // Add a green circle into the center
    cairo_new_path(cr);
    cairo_set_line_width(cr, 2.0f);
    cairo_set_source_rgb(cr, 0.0f, 1.0f, 0.0f);
    cairo_arc(cr, tag_x, tag_y, expectation, 0.0f, 2.0f * (float) M_PI);
    cairo_stroke(cr);

    // Add a magenta circle to the center of mass
    cairo_new_path(cr);
    cairo_set_line_width(cr, 2.0f);
    cairo_set_source_rgb(cr, 1.0f, 0.0f, 1.0f);
    cairo_arc(cr, center_x, center_y, 3.0f, 0.0f, 2.0f * (float) M_PI);
    cairo_stroke(cr);
    
    cairo_surface_flush(surface);
    cairo_destroy(cr);

    // Finally, add the anchors
    const double fn_size = (double) ((width < height) ? width : height) / 50.0;
    ls2_draw_anchors_to_cairo(surface, anchors, no_anchors, 0, fn_size);
}


extern void
ls2_cairo_write_pdf_inverted(const char* filename,
                             const float tag_x, const float tag_y,
			     const vector2 *anchors, const size_t no_anchors, 
			     const double *restrict result,
			     const uint16_t width, const uint16_t height,
			     const double center_x, const double center_y,
                             const double expectation, const double clip)
{
    cairo_surface_t *surface = cairo_pdf_surface_create(filename,
                                                        width, height);
    ls2_inverted_to_cairo_surface(surface, tag_x, tag_y, anchors, no_anchors,
                                  result, width, height, center_x, center_y,
                                  expectation, clip);
    cairo_surface_destroy(surface);
}


extern void __attribute__((__nonnull__))
ls2_cairo_write_png_inverted(const char* filename,
                             const float tag_x, const float tag_y,
			     const vector2 *anchors, const size_t no_anchors, 
			     const double *restrict result,
			     const uint16_t width, const uint16_t height,
			     const double center_x, const double center_y,
                             const double expectation, const double clip)
{
    cairo_surface_t *surface =
        cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
    ls2_inverted_to_cairo_surface(surface, tag_x, tag_y,
                                  anchors, no_anchors, 
				  result, width, height, center_x, center_y,
                                  expectation, clip);
    cairo_surface_write_to_png(surface, filename);
    cairo_surface_destroy(surface);
}
