# This file is part of LS² - the Localization Simulation Engine of FU Berlin.
#
# Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
# Stefan Adler, Malte Rohde, Jonathan Gunthermann
#
# LS² is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS² is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS².  If not, see <http://www.gnu.org/licenses/>.
#

EXTRA_DIST = algorithm/aml_algorithm.c \
	algorithm/aml_algorithm.h \
	algorithm/centroid_algorithm.c \
	algorithm/centroid_algorithm.h \
	algorithm/const_algorithm.c \
	algorithm/const_algorithm.h \
	algorithm/eminmax_w2_algorithm.c \
	algorithm/eminmax_w2_algorithm.h \
	algorithm/eminmax_w4_algorithm.c \
	algorithm/eminmax_w4_algorithm.h \
	algorithm/geo3_algorithm.c \
	algorithm/geo3_algorithm.h \
	algorithm/geon_algorithm.c \
	algorithm/geon_algorithm.h \
	algorithm/llsq_algorithm.c \
	algorithm/llsq_algorithm.h \
	algorithm/md_minmax_abs_algorithm.c \
	algorithm/md_minmax_abs_algorithm.h \
	algorithm/minmax_algorithm.c \
	algorithm/minmax_algorithm.h \
	algorithm/min_res1_bf_algorithm.c \
	algorithm/min_res1_bf_algorithm.h \
	algorithm/min_res2_bf_algorithm.c \
	algorithm/min_res2_bf_algorithm.h \
	algorithm/mle_gamma_algorithm.c \
	algorithm/mle_gamma_algorithm.h \
	algorithm/nllsq_algorithm.c \
	algorithm/nllsq_algorithm.h \
	algorithm/trilateration_algorithm.c \
	algorithm/trilateration_algorithm.h \
	algorithm/vble_algorithm.c \
	algorithm/vble_algorithm.h \
	algorithm/vble_opt_algorithm.c \
	algorithm/vble_opt_algorithm.h \
	algorithm/weighted_minmax_algorithm.c \
	algorithm/weighted_minmax_algorithm.h \
	algorithm/res_bruteforce_algorithm.c \
	algorithm/res_bruteforce_algorithm.h \
	backend/bits/circle.c \
	backend/bits/rectangle.c \
	estimator/crlb_malaney_algorithm.c \
	estimator/crlb_malaney_algorithm.h \
	estimator/crlb_qi_algorithm.c \
	estimator/crlb_qi_algorithm.h \
	estimator/crlb_so_algorithm.c \
	estimator/crlb_so_algorithm.h \
	estimator/crlb_zhao_algorithm.h \
	estimator/crlb_zhao_algorithm.c \
	estimator/gdop_algorithm.c \
	estimator/gdop_algorithm.h \
	error_model/bahillo_em.c error_model/bahillo_em.h \
	error_model/const_em.c error_model/const_em.h \
	error_model/eq_noise_em.c error_model/eq_noise_em.h \
	error_model/erlang_noise_em.c error_model/erlang_noise_em.h \
	error_model/gamma_noise_em.c error_model/gamma_noise_em.h \
	error_model/nd_noise_em.c error_model/nd_noise_em.h \
	error_model/nlosp_em.c error_model/nlosp_em.h \
	error_model/ray_noise_em.c error_model/ray_noise_em.h \
	util/util_circle.c \
	util/util_colors.c \
	util/util_math.c \
	util/util_matrix.c \
	util/util_median.c \
	util/util_misc.c \
	util/util_points.c \
	util/util_points_opt.c \
	util/util_random.c \
	util/util_triangle.c \
	util/util_vcircle.c \
	util/util_vector.c \
	avx_mathfun.h \
	sse_mathfun.h


dist_noinst_SCRIPTS = generate.py

DISTCLEANFILES = library.c ls2/library.h

BUILT_SOURCES = ls2/ls2-config.h library.c ls2/library.h

bin_PROGRAMS = ls2-run ls2-bounds ls2-diff ls2-h5-image
lib_LTLIBRARIES = libls2.la libls2be.la
nobase_include_HEADERS = ls2/ls2-config.h ls2/backend.h ls2/ls2.h \
	ls2/library.h ls2/output-variants.h

check_PROGRAMS = test-em

AM_CFLAGS = -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-align \
	-Wconversion -Wstrict-prototypes -Wmissing-prototypes -Wpacked \
	-Winline \
	-Werror -fdiagnostics-show-option \
	@ARCH_CFLAGS@ @RDRND_FLAGS@ -pthread \
	-ffast-math -fpredictive-commoning -ftree-vectorize

ls2_run_SOURCES  = shooter_main.c
ls2_run_CPPFLAGS =
ls2_run_LDADD    = $(LIBOBJS) libls2.la libls2be.la

ls2_bounds_SOURCES  = shooter_main.c
ls2_bounds_CPPFLAGS = -DESTIMATOR
ls2_bounds_LDADD    = $(LIBOBJS) libls2.la libls2be.la

ls2_diff_SOURCES  = ls2-diff.c
ls2_diff_CPPFLAGS =
ls2_diff_LDADD    = libls2be.la

ls2_h5_image_SOURCES  = ls2-h5-image.c
ls2_h5_image_CPPFLAGS =
ls2_h5_image_LDADD    = libls2be.la

libls2be_la_SOURCES = backend/backends.c backend/be_cairo.c \
		      backend/be_openexr.c backend/be_hdf5.c \
		      util/util_statistics.c ls2/util.h \
		      ls2/backend.h ls2/output-variants.h
libls2be_la_CPPFLAGS = $(CAIRO_CFLAGS) $(OPENEXR_CFLAGS)
libls2be_la_LDFLAGS  = -version $(LS2_LIBRARY_VERSION)
libls2be_la_LIBADD   = $(CAIRO_LIBS) $(OPENEXR_LIBS) $(LIBOBJS)

libls2_la_SOURCES  = shooter_run.c vector_shooter.h
libls2_la_CPPFLAGS = $(GSL_CFLAGS)
libls2_la_LDFLAGS  = -version $(LS2_LIBRARY_VERSION)
libls2_la_LIBADD   = $(GSL_LIBS) $(LIBOBJS)

test_em_SOURCES  = test-em.c vector_shooter.h ls2/ls2.h
test_em_CPPFLAGS =
test_em_LDADD    = $(LIBOBJS)

library.c: ${srcdir}/generate.py ${srcdir}/Makefile.am
	test -d ls2 || mkdir ls2
	${srcdir}/generate.py ${srcdir} $@
