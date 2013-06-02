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


 * SVN revision information:
 * @version $Revision$:
 * @author  $Author$:
 * @date    $Date$:
 * @URL     $URL$:
 */

#ifndef INCLUDED_LS2_BACKEND_H
#define INCLUDED_LS2_BACKEND_H

typedef enum ls2_output_format_t {
    OUTPUT_PNG,
    OUTPUT_PDF,
    OUTPUT_OPENEXR,
    NUM_OUTPUT_FORMATS,
} ls2_output_format_t;


enum ls2_output_format_t __attribute__((__const__))
get_output_format(const char *format);


extern void
ls2_write_locbased(ls2_output_format_t format, const char *filename,
		   const vector2 *anchors, const size_t num_anchors,
		   const float *results, const uint16_t width,
		   const uint16_t height);


extern void
ls2_write_inverted(ls2_output_format_t format, const char *filename,
		   const uint64_t runs, const float tag_x, float tag_y,
		   const vector2 *restrict anchors, const size_t num_anchors,
		   const uint64_t *restrict results, const uint16_t width,
		   const uint16_t height,
    		   const float center_x, float center_y);

/*
 *
 * BACKEND Functions.
 *
 */
extern void
ls2_cairo_write_pdf_locbased(const char* filename,
                             const vector2 *anchors, const size_t no_anchors,
                             const float* result, const uint16_t width,
		             const uint16_t height);

extern void
ls2_cairo_write_png_locbased(const char* filename,
                             const vector2 *anchors, const size_t no_anchors,
                             const float* result, const uint16_t width,
		             const uint16_t height);

extern void
ls2_cairo_write_pdf_inverted(const char* filename,
                             const float tag_x, const float tag_y,
			     const vector2 *anchors, const size_t no_anchors, 
			     const double *restrict result,
			     const uint16_t width, const uint16_t height,
			     const double center_x, const double center_y);

extern void
ls2_cairo_write_png_inverted(const char* filename,
                             const float tag_x, const float tag_y,
			     const vector2 *anchors, const size_t no_anchors, 
			     const double *restrict result,
			     const uint16_t width, const uint16_t height,
			     const double center_x, const double center_y);

extern void
ls2_cairo_write_png_phase_portrait(const char* filename,
				   const vector2 *anchors,
				   const size_t no_anchors, 
				   const float* dx, const float* dy,
				   const uint16_t width, const uint16_t height,
				   const uint16_t stride);

extern void
ls2_cairo_write_pdf_phase_portrait(const char* filename,
				   const vector2 *anchors,
				   const size_t no_anchors, 
				   const float* dx, const float* dy,
				   const uint16_t width, const uint16_t height,
				   const uint16_t stride);

extern void
ls2_cairo_write_pdf_diff(const char* filename, const vector2 *anchors,
                         const size_t no_anchors, const float* result,
                         const uint16_t width, const uint16_t height,
                         const double similar, const double dynamic);

extern void
ls2_cairo_write_png_diff(const char* filename, const vector2 *anchors,
                         const size_t no_anchors, const float* result,
                         const uint16_t width, const uint16_t height,
                         const double similar, const double dynamic);

extern void
ls2_openexr_write_locbased(const char* filename,
                           const vector2 *anchors, const size_t no_anchors,
                           const float* restrict result, const uint16_t width,
			   const uint16_t height);

extern void 
ls2_openexr_write_inverted(const char* filename,
                           const float tag_x, const float tag_y,
			   const vector2 *anchors, const size_t no_anchors,
                           const double *restrict result,
                           const uint16_t dim_x, const uint16_t dim_y,
    			   const double center_x, const double center_y);

extern int
ls2_openexr_write_diff(const char *filename, const vector2 *anchors,
                       const size_t no_anchors, const float *results,
                       const uint16_t width, const uint16_t height);

extern void
ls2_hdf5_write_locbased(const char *filename, const vector2 *anchors,
                        const size_t no_anchors, float **results,
                        const uint16_t width, const uint16_t height);

extern void 
ls2_hdf5_write_inverted(const char* filename,
                        const float tag_x, const float tag_y,
			const vector2 *restrict anchors, const size_t no_anchors,
                        const uint64_t *restrict result,
                        const uint16_t width, const uint16_t height,
    			const double center_x, const double center_y);

extern int
ls2_hdf5_write_diff(const char *filename, const vector2 *anchors,
                    const size_t no_anchors, const float *results[],
                    const uint16_t width, const uint16_t height);

extern int
ls2_hdf5_read_locbased(const char *filename, ls2_output_variant variant,
                       vector2 **anchors, size_t *no_anchors,
                       float **results, uint16_t *width, uint16_t *height);

extern int __attribute__((__nonnull__))
ls2_hdf5_read_inverted(const char *filename, float *tag_x, float *tag_y,
                       vector2 **anchors, size_t *no_anchors,
                       uint64_t **results, uint16_t *width, uint16_t *height,
                       double *center_x, double *center_y);

#endif
