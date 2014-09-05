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

EXTRA_DIST	= ls2.spec ls2.spec.in examples
SUBDIRS		= src tests
ACLOCAL_AMFLAGS	= -I m4

# These may be in the builddir too
BUILT_EXTRA_DIST = ChangeLog ls2.spec

${srcdir}/ls2.spec: ${srcdir}/ls2.spec.in ${srcdir}/configure.ac
	sed -e 's,@VERSION\@,${VERSION},g' ${srcdir}/ls2.spec.in > /tmp/ls2.spec ; \
	if test -f ${srcdir}/ls2.spec ; \
	then \
	    if cmp -s ${srcdir}/ls2.spec /tmp/ls2.spec ; \
	    then echo ls2.spec is unchanged ; rm /tmp/ls2.spec ; \
	    else rm -f ${srcdir}/ls2.spec ; mv /tmp/ls2.spec ${srcdir}/ls2.spec ; \
	    fi ; \
	else \
	    mv /tmp/ls2.spec ${srcdir}/ls2.spec ; \
	fi

ChangeLog:
	$(AM_V_GEN) if test -d "$(srcdir)/.git"; then \
	  (GIT_DIR=$(top_srcdir)/.git ./missing --run git log GLIB_2_20_0^^.. --stat) | fmt --split-only > $@.tmp \
	  && mv -f $@.tmp $@ \
	  || ($(RM) $@.tmp; \
	      echo Failed to generate ChangeLog, your ChangeLog may be outdated >&2; \
	      (test -f $@ || echo git-log is required to generate this file >> $@)); \
	else \
	  test -f $@ || \
	  (echo A git checkout and git-log is required to generate ChangeLog >&2 && \
	  echo A git checkout and git-log is required to generate this file >> $@); \
	fi


pkgconfigdir = $(datadir)/pkgconfig
pkgconfig_DATA = ls2.pc ls2be.pc

$(pkgconfig_DATA): config.status

html-local: Doxyfile
	test -d docs/api || mkdir -p docs/api
	doxygen Doxyfile

.PHONY: help
help:
	@echo
	@echo "Read INSTALL for installation instructions. To change the architecture, you"
	@echo "can reconfigure with appropriate settings or build with"
	@echo "  make ARCH_CFLAGS=\'options\´ clean all"
