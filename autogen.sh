#! /bin/sh
#
# Run this to generate all the initial scripts and makefiles, etc.
#
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

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

test -f $srcdir/src/generate.py || {
    echo "Error: Directors ${srcdir} does not look like the"
    echo "top-level source directory"
    exit 1
}

which libtool > /dev/null || {
    echo "You need to install libtool"
    exit 1
}

which autoreconf > /dev/null || {
    echo "You need to install automake and autoconf"
    exit 1
}

test -d $srcdir/m4 || mkdir $srcdir/m4
(cd $srcdir && autoreconf -i -s)
