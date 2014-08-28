#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

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

AC_PREREQ([2.67])
AC_INIT([LS2], [0.2], [heiko.will@fu-berlin.de])
AC_CONFIG_SRCDIR([src/vector_shooter.h])
AM_INIT_AUTOMAKE
AM_SILENT_RULES([yes])
AC_CONFIG_HEADERS([src/ls2/ls2-config.h])
AC_CONFIG_MACRO_DIR([m4])

# Initialize libtool
LT_PREREQ([2.2.6])
LT_INIT([disable-static,pic-only])
LS2_LIBRARY_VERSION=0:0:0
AC_SUBST([LS2_LIBRARY_VERSION])

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AM_PROG_CC_C_O
AC_PROG_MAKE_SET


# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AX_OPENCL(C)
PKG_CHECK_MODULES([GLIB2], [glib-2.0 >= 2.32])
PKG_CHECK_MODULES([GSL], [gsl >= 1.13])
PKG_CHECK_MODULES([CAIRO], [cairo >= 1.8])
PKG_CHECK_MODULES([OPENEXR], [OpenEXR >= 1.6])

AC_CHECK_LIB([hdf5], [H5Fcreate], [], [AC_MSG_FAILURE([You need to install libhdf5.so and hdf5.h])])
AC_CHECK_HEADERS([hdf5.h])

AX_OPENCL(C)

# Checks for header files.
# from AC_HEADER_ASSERT, changed for default no
AC_MSG_CHECKING([whether to enable assertions])
AC_ARG_ENABLE([assert],
  [AS_HELP_STRING([--enable-assert], [turn on assertions])],
  [ac_enable_assert=$enableval
   AS_IF(dnl
    [test "x$enableval" = xno],
      [AC_DEFINE([G_DISABLE_ASSERT], [1],
        [Define to 1 if assertions should be disabled.])],
    [test "x$enableval" != xyes],
      [AC_MSG_WARN([invalid argument supplied to --enable-assert])
      ac_enable_assert=yes])],
  [ac_enable_assert=no
   AC_DEFINE([G_DISABLE_ASSERT], [1])])
AC_MSG_RESULT([$ac_enable_assert])

AC_HEADER_STDBOOL
AC_HEADER_TIME
AC_CHECK_HEADERS([float.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T


dnl Check whether the compiler supports intrinsics.
AC_CHECK_HEADERS([immintrin.h], [],
    [AC_MSG_FAILURE([Does your compiler support intel''s intrinsics?])]
)

AC_MSG_CHECKING([whether CC supports -msse4.1])
my_save_cflags="$CFLAGS"
CFLAGS="$CFLAGS -msse4.1"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]],[[]])],
    [AC_MSG_RESULT([yes])]
    [cc_supports_sse41=yes],
    [AC_MSG_RESULT([no])]
    [cc_supports_sse41=no])

AC_MSG_CHECKING([for vector data type __m128])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#ifdef HAVE_IMMINTRIN_H
#include <immintrin.h>
#endif
]],[[
__m128 m = { 0, 0, 0, 0 };
m = _mm_round_ps(m, _MM_FROUND_TO_POS_INF);
]])],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])
AC_MSG_FAILURE([compiler does not support __m128])]
)
CFLAGS="$my_save_cflags"
unset my_save_cflags

AC_MSG_CHECKING([whether __m128 can be subscripted like an array])
my_save_cflags="$CFLAGS"
CFLAGS="$CFLAGS -msse4.1"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#ifdef HAVE_IMMINTRIN_H
#include <immintrin.h>
#endif
]],[[
__m128 m = { 0, 0, 0, 0 };
float x = m[0];
m[1] = 1.0f;
]])],
[AC_MSG_RESULT([yes])],
[AC_MSG_RESULT([no])
AC_MSG_FAILURE([CC cannot subscript __m128 like an array. Retry with an up-to-date compiler.])]
)
CFLAGS="$my_save_cflags"
unset my_save_cflags

my_save_cflags="$CFLAGS"
CFLAGS="$CFLAGS -mavx"
AC_MSG_CHECKING([whether CC supports -mavx])
AC_COMPILE_IFELSE([
     AC_LANG_PROGRAM([[#include <immintrin.h>]],
                     [[__m256 v = { 0, 0, 0, 0, 0, 0, 0, 0, };
                      _mm256_zeroupper();
                     ]])],
    [AC_MSG_RESULT([yes])]
    [cc_supports_avx=yes],
    [AC_MSG_RESULT([no])]
    [cc_supports_avx=no])
CFLAGS="$my_save_cflags"
unset my_save_cflags

if test "$cc_supports_avx" = "yes"
then
    AC_MSG_CHECKING([whether __m256 can be subscripted like an array])
    my_save_cflags="$CFLAGS"
    CFLAGS="-mavx"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <immintrin.h>]],[[
    __m256 m = { 0, 0, 0, 0, 0, 0, 0, 0 };
    float x = m[0];
    m[1] = 1.0f;
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
    AC_MSG_FAILURE([CC cannot subscript __m256 like an array. Retry with an up-to-date compiler.])]
)
    CFLAGS="$my_save_cflags"
    unset my_save_cflags
fi

my_save_cflags="$CFLAGS"
CFLAGS="$CFLAGS -mavx -mvzeroupper"
AC_MSG_CHECKING([whether CC supports -mvzeroupper])
AC_COMPILE_IFELSE([
     AC_LANG_PROGRAM([[#include <immintrin.h>]],
                     [[__m256 v = { 0, 0, 0, 0, 0, 0, 0, 0, };
                      _mm256_zeroupper();
                     ]])],
    [AC_MSG_RESULT([yes])]
    [cc_supports_vzeroupper=yes],
    [AC_MSG_RESULT([no])]
    [cc_supports_vzeroupper=no])
CFLAGS="$my_save_cflags"
unset my_save_cflags

my_save_cflags="$CFLAGS"
CFLAGS="$CFLAGS -mavx2"
AC_MSG_CHECKING([whether CC supports -mavx2])
AC_COMPILE_IFELSE([
     AC_LANG_PROGRAM([[#include <immintrin.h>]],
                     [[__m256 v = { 0, 0, 0, 0, 0, 0, 0, 0, };
                     ]])],
    [AC_MSG_RESULT([yes])]
    [cc_supports_avx2=yes],
    [AC_MSG_RESULT([no])]
    [cc_supports_avx2=no])
CFLAGS="$my_save_cflags"
unset my_save_cflags

if test "$cc_supports_avx2" = "yes"
then
    CFLAGS="$CFLAGS -mno-avx2"
fi

my_save_cflags="$CFLAGS"
CFLAGS="-mrdrnd"
AC_MSG_CHECKING([whether CC supports -mrdrnd])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]],[[]])],
    [AC_MSG_RESULT([yes])]
    [cc_supports_rdrnd=yes],
    [AC_MSG_RESULT([no])]
    [cc_supports_rdrnd=no])
CFLAGS="$my_save_cflags"
unset my_save_cflags

AC_ARG_WITH(rdrnd,
            [AS_HELP_STRING([--with-rdrnd],
                            [Use hardware generated random numbers. @<:@default=no@:>@])],
            [with_rdrnd=yes],
            [with_rdrnd=no])

if test -z "$RDRND_FLAGS" -a "x$with_rdrnd" = "xyes"
then
    RDRND_TESTS='rdrnd${EXEEXT}'
    if test "$cc_supports_rdrnd" = "yes"
    then
        RDRND_FLAGS="-mrdrnd"
    else
        AC_MSG_WARN([Use of rdrnd was requested but CC does not support -mrdrnd])
        RDRND_FLAGS=""
    fi
else
    RDRND_TESTS=""
    if test "$cc_supports_rdrnd" = "yes"
    then
        RDRND_FLAGS="-mno-rdrnd"
    else
        RDRND_FLAGS=""
    fi
fi
AC_SUBST(RDRND_FLAGS)
AC_SUBST(RDRND_TESTS)
AC_MSG_NOTICE([Building the library with rdrnd support?])
AC_MSG_NOTICE(AS_HELP_STRING([with rdrnd], [$with_rdrnd]))
AC_MSG_NOTICE(AS_HELP_STRING([rdrnd flag], [$RDRND_FLAGS]))

AC_ARG_WITH(arch,
            [AS_HELP_STRING([--with-arch],
                            [Build for architecture. @<:@default=native@:>@])],
            [],
            [with_arch=native])

if test -z "$ARCH_CFLAGS"
then
    ARCH_CFLAGS="-march=$with_arch"
fi
AC_SUBST(ARCH_CFLAGS)

# Checks for library functions.
AC_CHECK_FUNCS([memset sqrt strtoul])

# Test suite.
AC_CONFIG_TESTDIR([tests], [src])

AC_CONFIG_FILES([ls2be.pc
                 ls2.pc
                 Makefile
                 src/Makefile
		 tests/atlocal
                 tests/Makefile
                 Doxyfile])
AC_OUTPUT
