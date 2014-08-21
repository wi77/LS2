#  This file is part of LS² - the Localization Simulation Engine of FU Berlin.
#
#  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
#  Stefan Adler, Malte Rohde, Jonathan Gunthermann
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

R = R --vanilla
EXTRA_DIST = analyze.r

CLEANFILES = 
DISTCLEANFILES = 

BUILT_SOURCES = 

AM_CFLAGS = @ARCH_CFLAGS@ @RDRND_FLAGS@ -pthread -ffast-math -fpredictive-commoning -ftree-vectorize

check_PROGRAMS = $(RDRND_TEST) test-minres-bf test-algorithm test-error-model
EXTRA_PROGRAMS = rdrand

rdrand_SOURCES = rdrand.c
rdrand_CFLAGS = $(AM_CFLAGS)
rdrand_LDADD =

test_minres_bf_SOURCES = test-minres-bf.c
test_minres_bf_CPPFLAGS = -I${top_srcdir}/src -I../src
test_minres_bf_CFLAGS = @ARCH_CFLAGS@ -pthread -ffast-math -fpredictive-commoning -ftree-vectorize $(AM_CFLAGS) $(GLIB2_CFLAGS)
test_minres_bf_LDADD = $(GLIB2_LIBS) $(LIBOBJS)

test_error_model_SOURCES  = test-error-model.c \
	${top_srcdir}/src/vector_shooter.h \
	${top_srcdir}/src/ls2/ls2.h
test_error_model_CPPFLAGS = -I${top_srcdir}/src -I../src
test_error_model_CFLAGS   = -pthread $(AM_CFLAGS) $(GLIB2_CFLAGS)
test_error_model_LDADD    = $(GSL_LIBS) $(GLIB2_LIBS) $(LIBOBJS)

test_algorithm_SOURCES  = test-algorithm.c  \
	${top_srcdir}/src/vector_shooter.h \
	${top_srcdir}/src/ls2/ls2.h
test_algorithm_CPPFLAGS = -I${top_srcdir}/src -I../src
test_algorithm_CFLAGS   = -pthread $(AM_CFLAGS) $(GLIB2_CFLAGS)
test_algorithm_LDADD    = $(GSL_LIBS) $(GLIB2_LIBS) $(LIBOBJS)

