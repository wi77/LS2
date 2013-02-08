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

#ifndef INCLUDED_LS2_BE_TGA_H
#define INCLUDED_LS2_BE_TGA_H

#if defined(__GNUC_GNU_INLINE__) || defined(__GNUC_STD_INLINE__)

struct ls2_tga_image {
    uint8_t id;
    uint8_t color_map_type;
    uint8_t image_type;
    uint16_t cm_first_entry;
    uint16_t cm_length;
    uint8_t cm_bpp;
    uint16_t im_xorig;
    uint16_t im_yorig;
    uint16_t im_width;
    uint16_t im_height;
    uint8_t im_bpp;
    uint8_t im_flags;
    char image_id[20];
    uint8_t color_map[768];
    uint8_t image[]; //!< Size should be im_width * im_height.
} __attribute__((__packed__));

extern inline void __attribute__((__nonnull__,__leaf__,__gnu_inline__,__always_inline__,__artificial__))
ls2_tga_draw_pixel(struct ls2_tga_image *image, int color, int x, int y)
{
     if (0 <= x && x < image->im_width &&
	 0 <= y && y < image->im_height &&
         0 <= color && color <= 255) {
	  image->image[x + image->im_width * y] = (uint8_t) color;
     }
}

#else

struct ls2_tga_image;

extern void
ls2_tga_draw_pixel(struct ls2_tga_image *image, int color, int x, int y);

#endif

int __attribute__((__nonnull__))
ls2_tga_save(const char *filename, const struct ls2_tga_image *image);

void
ls2_tga_prepare_image(struct ls2_tga_image **image,
                      uint16_t width, uint16_t height);


void __attribute__((__nonnull__))
ls2_tga_free_image(struct ls2_tga_image *image);

void __attribute__((__nonnull__))
ls2_tga_set_color(struct ls2_tga_image * image,
		  uint8_t color, uint8_t r, uint8_t g, uint8_t b, uint8_t a);

extern void __attribute__((__nonnull__))
ls2_tga_write_locbased(const char* filename,
		       const size_t no_anchors, const vector2 *anchors,
		       const uint16_t width, const uint16_t height,
		       float* result);
#endif
