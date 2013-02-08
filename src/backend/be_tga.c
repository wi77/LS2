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
 * @author  $Author$:b
 * @date    $Date$:
 * @URL     $URL$:
 */
#include <immintrin.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "vector_shooter.h"

#include "backend/be_tga.h"


#if !defined(__GNUC_GNU_INLINE__) && !defined(__GNUC_STD_INLINE__)

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
    uint8_t image[];
} __attribute__((__packed__));

#endif

inline void
__attribute__((__nonnull__,__always_inline__,__artificial__))
ls2_tga_draw_pixel(struct ls2_tga_image *image, int x, int y, int color)
{
     if (0 <= x && x < image->im_width &&
         0 <= y && y < image->im_height &&
         0 <= color && color <= 255) {
	  image->image[x + image->im_width * y] = (uint8_t) color;
     }
}


#define DRAW_PIXEL(x, y) \
  ls2_tga_draw_pixel(image, x, y, color);

#include "bits/circle.c"
#include "bits/rectangle.c"

#if defined(__GNUC__) && (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 4)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

static inline
DECLARE_DRAW_CIRCLE(ls2_tga_draw_circle, struct ls2_tga_image *, int color)

static inline
DECLARE_DRAW_RECTANGLE(ls2_tga_draw_rectangle, struct ls2_tga_image *, int color)

#if defined(__GNUC__) && (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 4)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#endif


static void
ls2_tga_draw_result(struct ls2_tga_image *image,
		    uint16_t width, uint16_t height,
		    float* result)
{
    for (uint16_t y = 0; y < height; y++) {
        for (uint16_t x = 0; x < width; x++) {
            const int pos = x + width * y;
            const float sample = result[pos];
	    const int color = (sample < 250.0) ? (int) sample : 250;
	    ls2_tga_draw_pixel(image, x, y, color);
        }
    }
}


static void
ls2_tga_draw_anchors(struct ls2_tga_image *image,
		     const size_t no_anchors, const vector2 *anchors,
		     uint16_t offset)
{
    for (size_t i = 0; i < no_anchors; i++) {
	ls2_tga_draw_rectangle(image, image->im_height, image->im_width,
			       (uint16_t) (anchors[i].x + offset),
			       (uint16_t) anchors[i].y, 2, 2, 255);
    }
}


static struct ls2_tga_image * __attribute__((__malloc__))
allocate_and_setup_header(uint16_t width, uint16_t height)
{
    
    const size_t sz = (size_t)(width * height) * sizeof(uint8_t);
    struct ls2_tga_image *image = malloc(sizeof(struct ls2_tga_image) + sz);
    if (image != NULL) {
	/* Fill the header */
	image->id = 20;
	image->color_map_type = 1;
	image->image_type = 1;
	image->cm_first_entry = 0;
	image->cm_length = 256;
	image->cm_bpp = 24;
	image->im_xorig = 0;
	image->im_yorig = 0;
	image->im_width = width;
	image->im_height = height;
	image->im_bpp = 8;
	image->im_flags = 32;
	strncpy(image->image_id, "created by ls2     ", 20);
    }
    return image;
}



void __attribute__((__nonnull__))
ls2_tga_set_color(struct ls2_tga_image * image,
		  uint8_t color, uint8_t r, uint8_t g, uint8_t b,
                  uint8_t a __attribute__((__unused__)))
{
    image->color_map[color * 3 + 2] = r;
    image->color_map[color * 3 + 1] = g;
    image->color_map[color * 3 + 0] = b;
    // ignore a.
}



void
ls2_tga_prepare_image(struct ls2_tga_image **image,
		      uint16_t width, uint16_t height)
{

    *image = allocate_and_setup_header(width, height);


    uint8_t i=0;    
    // build color table
    do {
        if (i <= (VERY_GOOD_VALUES_COLOR)) {
	    ls2_tga_set_color(*image, i, (uint8_t)(5 * i), 255, (uint8_t)(5 * i), 0);
        } else {
	    uint8_t luminance = (uint8_t) (250.0f - ((float) i - 50.0f) * (250.0f/200.0f));
	    ls2_tga_set_color(*image, i, luminance, luminance, luminance, 0);
        }
        i++;  
    } while (i<250);
   
    // Max failure in blue
    ls2_tga_set_color(*image, 250, 0, 0, 255, 0);

    // No hit in white.
    ls2_tga_set_color(*image, 251, 255, 255, 255, 0);

    // Max failure in blue
    ls2_tga_set_color(*image, 252, 0, 0, 255, 0);

    // Center of localisation in magenta.
    ls2_tga_set_color(*image, 253, 255, 0, 255, 0);

    // Tag location in yellow
    ls2_tga_set_color(*image, 254, 255, 255, 0, 0);

    // Anchors in red
    ls2_tga_set_color(*image, 255, 255, 0, 0, 0);
}



void __attribute__((__nonnull__))
ls2_tga_free_image(struct ls2_tga_image *image)
{
    free(image);
}



// saves an array of pixels as a TGA image
int __attribute__((__nonnull__))
ls2_tga_save(const char *filename, const struct ls2_tga_image *image)
{
    /* Write a TGA image */
    ssize_t written;
    int fd = open(filename, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
    if (fd == -1) goto error1;
    const size_t sz = sizeof(struct ls2_tga_image) +
        (size_t)(image->im_width * image->im_height) * sizeof(uint8_t);
    written = write(fd, image, sz);
    if (written != (ssize_t) sz) goto error;

    close(fd);
    return 0;

  error:
    close(fd);
  error1:
    perror("ls2_tga_save()");
    return -1;
}



extern void __attribute__((__nonnull__))
ls2_tga_write_locbased(const char* filename,
		       const size_t no_anchors, const vector2 *anchors,
		       const uint16_t width, const uint16_t height,
		       float* result)
{
    struct ls2_tga_image *image;
    ls2_tga_prepare_image(&image, width, height);
    if (image == NULL) {
        perror("ls2_tga_write_locbased(): prepare()");
        exit(EXIT_FAILURE);
    }

    ls2_tga_draw_result(image, width, height, result);
    ls2_tga_draw_anchors(image, no_anchors, anchors, 0);
    ls2_tga_save(filename, image);

  // finish:
    ls2_tga_free_image(image);
}


