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

#undef  _XOPEN_SOURCE
#define _XOPEN_SOURCE 600

#undef  _GNU_SOURCE
#define _GNU_SOURCE

#ifdef HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#include <stdlib.h>
#include <stdint.h>

/* By defining the macro H5_NO_DEPRECATED_SYMBOLS, we force the use of
 * the version 1.8 API regardless of the configuration of hdf5
 */
#define H5_NO_DEPRECATED_SYMBOLS 1
#include <hdf5.h>

#include "ls2/ls2.h"
#include "ls2/backend.h"


static const char *ls2_hdf5_variant[] = {
#undef  LS2OUT_VARIANT
#define LS2OUT_VARIANT(tag, name, h5name) h5name,
#include "ls2/output-variants.h"
    NULL
};



static const char *
ls2_hdf5_variant_name(ls2_output_variant variant)
{
    if (variant < NUM_VARIANTS)
        return ls2_hdf5_variant[variant];
    else
        return NULL;
}




static void
ls2_hdf_write_anchors(hid_t file_id, const vector2 *anchors, size_t no_anchors)
{
    hid_t dataset, dataspace;
    hsize_t dims[2] = { no_anchors, 2 };
    
    dataspace = H5Screate_simple(2, dims, NULL);
    dataset = H5Dcreate(file_id, "/Anchors", H5T_NATIVE_FLOAT,
                        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             anchors);
    H5Sclose(dataspace);
    H5Dclose(dataset);
}



void
ls2_hdf5_write_locbased(const char *filename, const vector2 *anchors,
                        const size_t no_anchors, float **results,
                        const uint16_t width, const uint16_t height)
{
    hid_t file_id, grp, dataset, dataspace, plist_id;
    hsize_t dims[2];
    hsize_t chunk_dims[2] = { width, height };

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    grp = H5Gcreate(file_id, "/Result", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    ls2_hdf_write_anchors(file_id, anchors, no_anchors);

    dims[0] = height;
    dims[1] = width;
    for (ls2_output_variant k = 0; k < NUM_VARIANTS; k++) {
        if (results[k] == NULL)
            continue;
        char name[256];
        dataspace = H5Screate_simple(2, dims, NULL);
	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist_id, 2, chunk_dims);
	H5Pset_deflate (plist_id, 9);
        snprintf(name, 256, "/Result/%s", ls2_hdf5_variant_name(k));
        dataset = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT,
                            dataspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 results[k]);
        H5Sclose(dataspace);
        H5Dclose(dataset);
    }
    H5Gclose(grp);
    H5Fclose(file_id);
}





static int __attribute__((__nonnull__))
ls2_hdf5_read_anchors(hid_t file, vector2** anchors, size_t *no_anchors)
{
    hid_t dataset, dataspace, memspace;
    hsize_t dims[2];
    int rank;

    dataset = H5Dopen2(file, "/Anchors", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2) {
        fprintf(stderr, "/Anchors wrong rank %d\n", (int) rank);
        return -1;
    }
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (dims[1] != 2) {
        fprintf(stderr, "/Anchors dimension (%d, %d)\n",
                (int) dims[0], (int) dims[1]);
        return -1;
    }
    *no_anchors = (size_t) dims[0];
    *anchors = (vector2 *) calloc(*no_anchors, sizeof(vector2));
    memspace = H5Screate_simple(2, dims, NULL);
    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT,
            *anchors);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);

    return 0;
}




int __attribute__((__nonnull__(1,5,6,7)))
ls2_hdf5_read_locbased(const char *filename, ls2_output_variant variant,
                       vector2 **anchors, size_t *no_anchors,
                       float **results, uint16_t *width, uint16_t *height)
{
    hid_t file, dataset, dataspace, memspace;
    hsize_t dims[2];
    int rank, ret;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Read the anchors.
    if (anchors != NULL && no_anchors != NULL) {
        ret = ls2_hdf5_read_anchors(file, anchors, no_anchors);
        if (ret < 0)
            return ret;
    }

    // Read the errors.
    char name[256];
    snprintf(name, 256, "/Result/%s", ls2_hdf5_variant_name(variant));
    dataset = H5Dopen2(file, name, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2) {
        fprintf(stderr, "%s wrong rank %d\n", name, (int) rank);
        return -1;
    }
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    *width = (uint16_t) dims[1];
    *height = (uint16_t) dims[0];
    *results = (float *) calloc((size_t) (*height * *width), sizeof(float));
    memspace = H5Screate_simple(2, dims, NULL);
    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT,
            *results);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);

    H5Fclose(file);
    return 0;
}




void
ls2_hdf5_write_inverted(const char *filename,
			const float tag_x, const float tag_y,
			const vector2 *restrict anchors, const size_t no_anchors,
			const uint64_t *restrict result,
			const uint16_t width, const uint16_t height,
			const double center_x, const double center_y)
{
    hid_t file_id, grp, dataset, dataspace, plist_id;
    hsize_t dims[2];
    hsize_t chunk_dims[2] = { width, height };

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    grp = H5Gcreate(file_id, "/Result", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    ls2_hdf_write_anchors(file_id, anchors, no_anchors);

    dims[0] = 2;
    dims[1] = 1;
    dataspace = H5Screate_simple(2, dims, NULL);
    dataset = H5Dcreate(file_id, "/Tag", H5T_NATIVE_FLOAT,
                        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    float tag[2] = {tag_x, tag_y};
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tag);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    dims[0] = 2;
    dims[1] = 1;
    dataspace = H5Screate_simple(2, dims, NULL);
    dataset = H5Dcreate(file_id, "/Result/Center", H5T_NATIVE_DOUBLE,
                        dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double center[2] = { center_x, center_y};
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, center);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    dims[0] = height;
    dims[1] = width;
    dataspace = H5Screate_simple(2, dims, NULL);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 2, chunk_dims);
    H5Pset_deflate (plist_id, 9);
    // BUG: should be a native type, but what is uint64_t?
    dataset = H5Dcreate(file_id, "/Result/Frequencies", H5T_STD_U64LE,
			dataspace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
	     result);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Gclose(grp);
    H5Fclose(file_id);
}





int __attribute__((__nonnull__))
ls2_hdf5_read_inverted(const char *filename, float *tag_x, float *tag_y,
                       vector2 **anchors, size_t *no_anchors,
                       uint64_t **results, uint16_t *width, uint16_t *height,
                       double *center_x, double *center_y)
{
    hid_t file, dataset, dataspace, memspace;
    hsize_t dims[2];
    int rank, ret;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // Read the anchors.
    if ((ret = ls2_hdf5_read_anchors(file, anchors, no_anchors)) < 0)
        return ret;

    // Read the tag position
    dataset = H5Dopen2(file, "/Tag", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2) {
        fprintf(stderr, "/Tag wrong rank %d\n", (int) rank);
        return -1;
    }
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    // TODO: Validate values.
    float tag[2];
    memspace = H5Screate_simple(2, dims, NULL);
    H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, tag);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    *tag_x = tag[0];
    *tag_y = tag[1];

    // Read the centroid of all estimates
    dataset = H5Dopen2(file, "/Result/Center", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2) {
        fprintf(stderr, "/Result/Center wrong rank %d\n", (int) rank);
        return -1;
    }
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    // TODO: Validate values.
    double center[2];
    memspace = H5Screate_simple(2, dims, NULL);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
            center);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    *center_x = center[0];
    *center_y = center[1];

    // Read the frequencies
    dataset = H5Dopen2(file, "/Result/Frequencies", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 2) {
        fprintf(stderr, "/Result/Frequencies wrong rank %d\n", (int) rank);
        return -1;
    }
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    *width = (uint16_t) dims[1];
    *height = (uint16_t) dims[0];
    *results = (uint64_t*) calloc((size_t) (*height * *width),
                                  sizeof(uint64_t));
    memspace = H5Screate_simple(2, dims, NULL);
    H5Dread(dataset, H5T_STD_U64LE, memspace, dataspace, H5P_DEFAULT,
            *results);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);

    H5Fclose(file);
    return 0;
}
