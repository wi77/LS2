#! /usr/bin/python
# -*- coding: utf8 -*-
#
# Automatically generate dispatch tables for the library
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

import re
import os
import sys
import string


def find_algorithms():
    """Find all algorithm files."""
    algfilter = re.compile("(.*)_algorithm.c$")
    dirlist =  os.listdir(sys.argv[1] + "/algorithm")
    algs = [ algfilter.match(x).groups()[0] for x in dirlist if algfilter.match(x) ]
    algs.sort()
    return algs


def find_error_models():
    """Find all error model files."""
    emfilter = re.compile("(.*)_em.c$")
    dirlist = os.listdir(sys.argv[1] + "/error_model")
    ems = [ emfilter.match(x).groups()[0] for x in dirlist if emfilter.match(x) ]
    ems.sort()
    return ems


def find_estimators():
    """Find all bound estimator implementations."""
    estfilter = re.compile("(.*)_algorithm.c$")
    dirlist = os.listdir(sys.argv[1] + "/estimator")
    ests = [ estfilter.match(x).groups()[0] for x in dirlist if estfilter.match(x) ]
    ests.sort()
    return ests


def algorithm_file(alg):
    return sys.argv[1] + "/algorithm/" + alg + "_algorithm.c"


def find_algorithm_name(alg):
    """Find the name of an algorithm in file alg_file."""
    expr = re.compile(".*@algorithm_name: (.*) \*/.*")
    for line in open(algorithm_file(alg)):
        m = expr.match(line)
        if m:
            return m.groups()[0]
    return alg


def estimator_file(est):
    return sys.argv[1] + "/estimator/" + est + "_algorithm.c"


def find_estimator_name(est):
    """Find the name of an estimator in file est_file."""
    expr = re.compile(".*@algorithm_name: (.*) \*/.*")
    for line in open(estimator_file(est)):
        m = expr.match(line)
        if m:
            return m.groups()[0]
    return est


def error_model_file(em):
    return sys.argv[1] + "/error_model/" + em + "_em.c"


def find_error_model_name(em):
    """Find the name of an error model in file em_file."""
    expr = re.compile(".*@error_model_name: (.*) \*/.*")
    for line in open(error_model_file(em)):
        m = expr.match(line)
        if m:
            return m.groups()[0]
    return em


def find_command_line_arguments(f):
    """Checks whether an algorithm, estimator, or error model defines command
    line parameters."""
    expr = re.compile(".*GOptionEntry (.*)\[\].*")
    for line in open(f):
        m = expr.match(line)
        if m:
            return m.groups()[0]
    return None

algs = find_algorithms()
ems = find_error_models()
ests = find_estimators()

lib = open("library.c", "w")
head = open("ls2/library.h", "w")

head.writelines([ '/* This file was automatically generated. Do not edit! */\n',
                 '\n',
                ])
lib.writelines([ '/* This file was automatically generated. Do not edit! */\n',
                 '\n'
               ])

lib.writelines([ '#include "algorithm/' + alg + '_algorithm.c"\n' for alg in algs ])
lib.writelines([ '\n' ])
lib.writelines([ '#include "error_model/' + em + '_em.c"\n' for em in ems ])

lib.writelines([ '\n' ])
lib.writelines([ '#include "estimator/' + est + '_algorithm.c"\n' for est in ests ])

head.writelines([ 'typedef enum algorithm_t {\n' ])
head.writelines([ '    ALG_' + alg.upper() + ',\n' for alg in algs ])
head.writelines([ "} algorithm_t;\n",
                  "\n"
                ])

head.write('\n#define ALGORITHMS \"' + string.join([ '\\\"' + alg.replace('_', '-') + '\\\"' for alg in algs], ', ') + '\"\n\n')
lib.writelines([ '\nconst char* const algorithm_short_name[] = {\n' ])
lib.writelines([ '    "' + alg.replace('_', '-') + '",\n' for alg in algs ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "algorithm_t get_algorithm_by_name(const char *name)\n",
                 "{\n",
                 "    int i = 0;\n",
                 "    while (algorithm_short_name[i] != NULL) {\n",
                 "        if (strcmp(name, algorithm_short_name[i]) == 0) return i;\n",
                 "        i++;\n",
                 "    }\n",
                 "    return -1;\n",
                 "}\n",
                 "\n" ])

lib.writelines([ 'const char* const algorithm_name[] = {\n' ])
lib.writelines([ '    "' + find_algorithm_name(alg) + '",\n' for alg in algs ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "const char * const * get_algorithms(void)\n",
                 "{\n",
                 "    return algorithm_name;\n",
                 "}\n",
                 "\n" ])

lib.writelines([ 'static inline void __attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))\n',
                 'algorithm(algorithm_t alg,\n',
                 '          const VECTOR *restrict vx, const VECTOR *restrict vy, const VECTOR *restrict r,\n',
                 '          size_t no_anchors, int width, int height,\n',
                 '          VECTOR *restrict resx, VECTOR *restrict resy)\n',
                 '{\n',
                 '    switch (alg) {\n',
                ])
lib.writelines([ '    case ALG_' + alg.upper() + ':\n        ' + alg + '_run(vx, vy, r, no_anchors, width, height, resx, resy);\n        break;\n' for alg in algs ])
lib.writelines([ "    }\n",
                 "}\n",
                 "\n",
               ])


head.writelines([ 'typedef enum estimator_t {\n' ])
head.writelines([ '    EST_' + est.upper() + ',\n' for est in ests ])
head.writelines([ "} estimator_t;\n",
                  "\n"
                ])

head.write('\n#define ESTIMATORS \"' + string.join([ '\\\"' + est.replace('_', '-') + '\\\"' for est in ests], ', ') + '\"\n\n')
lib.writelines([ '\nconst char* const estimator_short_name[] = {\n' ])
lib.writelines([ '    "' + est.replace('_', '-') + '",\n' for est in ests ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "estimator_t get_estimator_by_name(const char *name)\n",
                 "{\n",
                 "    int i = 0;\n",
                 "    while (estimator_short_name[i] != NULL) {\n",
                 "        if (strcmp(name, estimator_short_name[i]) == 0) return i;\n",
                 "        i++;\n",
                 "    }\n",
                 "    return -1;\n",
                 "}\n",
                 "\n" ])

lib.writelines([ 'const char* const estimator_name[] = {\n' ])
lib.writelines([ '    "' + find_estimator_name(est) + '",\n' for est in ests ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "const char * const * get_estimators(void)\n",
                 "{\n",
                 "    return estimator_name;\n",
                 "}\n",
                 "\n" ])

lib.writelines([ 'static inline float __attribute__((__always_inline__,__pure__))\n',
                 'estimate(estimator_t est, const vector2 *anchor, const size_t num_anchors, const vector2 *location)\n',
                 '{\n',
                 '    switch (est) {\n',
                ])
lib.writelines([ '    case EST_' + est.upper() + ':\n        return ' + est + '_run(anchor, num_anchors, location);\n        break;\n' for est in ests ])
lib.writelines([ "    }\n",
                 "    return NAN;\n",
                 "}\n",
                 "\n",
               ])


head.writelines([ 'typedef enum error_model_t {\n' ])
head.writelines([ '    EM_' + em.upper() + ',\n' for em in ems ])
head.writelines([ "} error_model_t;\n",
                  "\n"
                ])

head.write('\n#define ERROR_MODELS \"' + string.join([ '\\\"' + em.replace('_', '-') + '\\\"' for em in ems ], ', ') + '\"\n\n')
lib.writelines([ 'const char* const error_model_short_name[] = {\n' ])
lib.writelines([ '    "' + em.replace('_', '-') + '",\n' for em in ems ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "error_model_t get_error_model_by_name(const char *name)\n",
                 "{\n",
                 "    int i = 0;\n",
                 "    while (error_model_short_name[i] != NULL) {\n",
                 "        if (strcmp(name, error_model_short_name[i]) == 0) return i;\n",
                 "        i++;\n",
                 "    }\n",
                 "    return -1;\n",
                 "}\n",
                 "\n" ])

lib.write("const char* const error_model_name[] = {\n")
lib.writelines([ '    "' + find_error_model_name(em) + '",\n' for em in ems ])
lib.writelines([ "    NULL\n",
                 "};\n",
                 "\n",
                 "const char * const * get_error_models(void)\n",
                 "{\n",
                 "    return error_model_name;\n",
                 "}\n",
                 "\n" ])


lib.writelines([ 'static void\n',
                 'error_model_setup(error_model_t model, const vector2 *anchors, size_t nanchors)\n',
                 '{\n',
                 '    switch (model) {\n',
                ])
lib.writelines([ '    case EM_' + em.upper() + ':\n        ' + em + '_setup(anchors, nanchors);\n        break;\n' for em in ems ])
lib.writelines([ "    }\n",
                 "}\n",
                 "\n",
               ])

lib.writelines([ 'static inline void __attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__(2,3,9)))\n',
                 'error_model(error_model_t model, __m128i *restrict seed, const VECTOR *restrict dist,\n',
                 '            const VECTOR *restrict vx, const VECTOR *restrict vy, size_t no_anchors,\n'
                 '            const VECTOR tagx, const VECTOR tagy, VECTOR *restrict result)\n',
                 '{\n',
                 '    switch (model) {\n',
                ])
lib.writelines([ '    case EM_' + em.upper() + ':\n        ' + em + '_error(seed, no_anchors, dist, vx, vy, tagx, tagy, result);\n        break;\n' for em in ems ])
lib.writelines([ "    }\n",
                 "}\n",
                 "\n",
               ])

# Collect and write out the command line parameter tables if each file
# defines one.
lib.write("extern void __attribute__((__nonnull__))\nls2_add_algorithm_option_groups(GOptionContext *context)\n{\n")
for alg in algs:
    args = find_command_line_arguments(algorithm_file(alg))
    if args:
        lib.write("    ls2_add_%s_option_group(context);\n" % alg)
lib.write("\n}\n\n")
lib.write("extern void __attribute__((__nonnull__))\nls2_add_error_model_option_groups(GOptionContext *context)\n{\n")
for em in ems:
    args = find_command_line_arguments(error_model_file(em))
    if args:
        lib.write("    ls2_add_%s_option_group(context);\n" % em)
lib.write("\n}\n\n")
lib.write("extern void __attribute__((__nonnull__))\nls2_add_estimator_option_groups(GOptionContext *context)\n{\n")
for est in ests:
    args = find_command_line_arguments(estimator_file(est))
    if args:
        lib.write("    ls2_add_%s_option_group(context);\n" % est)
lib.write("\n}\n")

head.flush()
head.close()
lib.flush()
lib.close()
