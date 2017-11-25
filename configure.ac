dnl              xa_stats : statistics
dnl

#AC_C_BIGENDIAN

DEBUG_CXX=g++
DEBUG_CC=gcc
NORMAL_CXX=g++
NORMAL_CC=gcc
# the compiler will be decided later when we see the +/- debug switch, for the time use this
CXX=${NORMAL_CXX}
CC=${NORMAL_CC}

#MY_CFLAGS='-finline-functions -funroll-all-loops -O3 -Wall -ansi -pedantic -W' # will be added later
#MY_CFLAGS='-march=core2 -mtune=generic -O3 -pipe -Wall -W' # will be added later
# the best is this with g++-fsf-4.6
MY_LDFLAGS=""
MY_CFLAGS='-Wall -W' # will be added later
#MY_CFLAGS='-O3 -Wall -ansi -pedantic -W' # will be added later
MY_CPPFLAGS=${MY_CFLAGS}

AC_INIT([xa_stats], [3.0.0], [andreashad2@gmail.com ex livantes@soi.city.ac.uk])
AC_CONFIG_SRCDIR([src/C/anova.c])
AM_INIT_AUTOMAKE([1.10 no-define])
AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_PREREQ([2.69]) # latest autoconf version

AC_CANONICAL_HOST
#AC_CANONICAL_TARGET

LT_INIT

m4_pattern_allow([AC_CONFIG_MACRO])
m4_pattern_allow([AS_AC_EXPAND])
m4_pattern_allow([AX_CC_MAXOPT])
m4_pattern_allow([AX_BOOST_BASE])
m4_pattern_allow([AC_MSG_ERROR])
m4_pattern_allow([AC_MSG_WARN])
m4_pattern_allow([AM_PROG_LIBTOOL])
m4_pattern_allow([AC_DEFINE])

echo "You have host: ${host}, cpu: ${host_cpu}, os: ${host_os}, vendor: ${host_vendor}"
echo "Target : ${target}, cpu: ${target_cpu}, os: ${target_os}, vendor: ${target_vendor}"

#AC_PROG_RANLIB
AC_PROG_LIBTOOL

#AM_GNU_GETTEXT([external])
#AM_GNU_GETTEXT_VERSION(0.15)
#AM_ICONV

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADER(tr1/unordered_map,[AC_DEFINE([HAVE_TR1],[],["Have tr1"])],[])
AC_CHECK_HEADER(unordered_map,[AC_DEFINE([HAVE_CXX0X],[],["Have C++0x"])],[])

AC_DEFUN([AC_ARG_DISABLE], [
    disablestring=`echo "$1" | sed 's/_/-/g'` 
    AC_MSG_CHECKING([--disable-$disablestring])
    AC_ARG_ENABLE($1, [$2], 
    [
        if test "$enableval" = yes; then
            ac_cv_use_$1='$3=yes'
        else
            ac_cv_use_$1='$3=no'
        fi
    ],
    [
        ac_cv_use_$1='$3='$DEFAULT_$3
    ])

    eval "$ac_cv_use_$1"

    if test "$$3" = "no"; then
        AC_MSG_RESULT(yes)
    else
        AC_MSG_RESULT(no)
    fi
])

# debug
# by default threads are on, so actions are reversed
#pthreads=yes
#AC_ARG_ENABLE(pthreads,[AC_HELP_STRING([--disable-threads],[disable multi-threading using POSIX threads (pthreads) which must be installed in the system already - this provides huge savings in running time and should be disabled only if you all else fails. By default pthreads are enabled.])], [pthreads=$enableval], [pthreads=yes])
#AM_CONDITIONAL(PTHREADS, test "$pthreads" = "yes")
# by default warn-all are on, so actions are reversed
AC_ARG_ENABLE(werror, [AC_HELP_STRING([--disable-warn_all],[DO NOT stop at each warning during compilation.])], [werror=$enableval], [werror=yes])
AM_CONDITIONAL(WERROR, test "$werror" = "yes")
AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[disable ALL optimisations and use -g flag for compilation (produces executable for to be inspected by gdb - unlikely to be used by regular user)])], [debug=$enableval], [debug=no])
AM_CONDITIONAL(DEBUG, test "$debug" = "yes")
AC_ARG_ENABLE(profiling, [AC_HELP_STRING([--enable-profiling],[disable ALL optimisations and use -pg flag for profiling the executable with gprof. Unlikely to be needed by regular user.])], [profiling=$enableval], [profiling=no])
AM_CONDITIONAL(PROFILING, test "$profiling" = "yes")
AC_ARG_ENABLE(verbose, [AC_HELP_STRING([--enable-verbose],[print a lot of info regarding insertion/deletion of nodes, etc. WARNING: if the input files are large dont turn this on!])], [verbose=$enableval], [verbose=no])
AM_CONDITIONAL(VERBOSE, test "$verbose" = "yes")
AC_ARG_ENABLE(valgrind, [AC_HELP_STRING([--enable-valgrind],[turn off all optimisations so that valgrind can do its job, don't forget to use the GLIB-friendly valgrind switches (export G_SLICE='always-malloc' and export G_DEBUG='gc-friendly,resident-modules', then 'valgrind --dsymutil=yes --tool=memcheck --leak-check=full', note that '--dsymutil=yes' is a special switch for idiotic MACOSX which removes all debug symbols.])], [valgrind=$enableval], [valgrind=no])
AM_CONDITIONAL(VALGRIND, test "$valgrind" = "yes")

# input data file consistency checks warnings
AC_ARG_ENABLE(data_file_consistency_checks, [AC_HELP_STRING([--enable-data_file_consistency_checks],[keep any optimisations but check input data file consistency, e.g. are all edges mentioned in the edges data file are also contained in the input edges file?)])], [data_file_consistency_checks=$enableval], [data_file_consistency_checks=no])
AM_CONDITIONAL(DATA_FILE_CONSISTENCY_CHECKS, test "$data_file_consistency_checks" = "yes")

# optimisations, these will not have an effect if the --enable-debug above is set
AC_ARG_ENABLE(sse, [AC_HELP_STRING([--disable-sse],[disable SSE optimizations])], have_sse=$enableval, have_sse=yes)
if test "$have_sse" = "yes"; then have_sse2="yes"; fi

AC_ARG_ENABLE(sse2, [AC_HELP_STRING([--disable-sse2],[disable SSE2 optimizations])], have_sse2=$enableval, have_sse2=yes)
AM_CONDITIONAL(HAVE_SSE2, test "$have_sse2" = "yes")

AC_ARG_ENABLE(altivec, [AC_HELP_STRING([--disable-altivec],[disable Altivec optimizations])], have_altivec=$disableval, have_altivec=yes)
AM_CONDITIONAL(HAVE_ALTIVEC, test "$have_altivec" = "yes")

AC_ARG_ENABLE(pipe, [AC_HELP_STRING([--disable-pipe],[disable pipe optimizations])], have_pipe=$disableval, have_pipe=yes)
AM_CONDITIONAL(HAVE_PIPE, test "$have_pipe" = "yes")

AC_ARG_ENABLE(ftree-vectorize, [AC_HELP_STRING([--disable-ftree-vectorize],[disable Tree Vectorize optimizations])], have_ftree_vectorize=$enableval, have_ftree_vectorize=yes)
AM_CONDITIONAL(HAVE_TREE_VECTORIZE, test "$have_ftree_vectorize" = "yes")

AC_ARG_ENABLE(avx, [AC_HELP_STRING([--enable-avx],[enable AVX optimizations])], have_avx=$enableval, have_avx=no)
AM_CONDITIONAL(HAVE_AVX, test "$have_avx" = "yes")

AC_ARG_ENABLE(neon, [AC_HELP_STRING([--enable-neon],[enable ARM NEON optimizations])], have_neon=$enableval, have_neon=no)
AM_CONDITIONAL(HAVE_NEON, test "$have_neon" = "yes")

HAVE_PTHREADS=""
if test "$pthreads" = "yes"; then
	AX_PTHREAD([
		HAVE_PTHREADS="-DHAVE_PTHREADS"
	],[
		AC_MSG_WARN(
			[failed to find pthreads support, will not use pthreads, the program will not be multi-threaded.]
		)
	])
fi

# don't do the optimisation flags if on debug mode
OPTIMISATION_DETECTED=""
DEBUG_FLAGS=""
if test "$debug" = "no" && test "$profiling" = "no"; then
	CXX=${NORMAL_CXX}
	CC=${NORMAL_CC}
else
	# debug
	CXX=${DEBUG_CXX}
	CC=${DEBUG_CC}
fi

# check for C++ preprocessor and compiler and the library compiler
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CPP
AC_PROG_CC
AM_PROG_CC_C_O

if test "$debug" = "no" && test "$profiling" = "no"; then
	AX_CC_MAXOPT
	case "${ax_cv_c_compiler_vendor}" in
	    gnu|intel)
		# SSE/SSE2
		if test "$have_sse2" = "yes" -a "x$SSE2_CFLAGS" = x; then
		    AX_CHECK_COMPILER_FLAGS(-msse2, [SSE2_CFLAGS="-msse2"],
			[AC_MSG_WARN([Need a version of gcc with -msse2])])
		fi

		# Tree vectorise
		if test "$have_ftree_vectorize" = "yes" -a "x$TREE_VECTORIZE_CFLAGS" = x; then
		    AX_CHECK_COMPILER_FLAGS(-ftree-vectorize, [TREE_VECTORIZE_CFLAGS="-ftree-vectorize"],
			[AC_MSG_WARN([Need a version of gcc with -ftree-vectorize])])
		fi

		# -pipe
		if test "$have_pipe" = "yes" -a "x$PIPE_CFLAGS" = x; then
		    AX_CHECK_COMPILER_FLAGS(-pipe, [PIPE_CFLAGS="-pipe"],
			[AC_MSG_WARN([Need a version of gcc with -pipe])])
		fi

		# AVX
		if test "$have_avx" = "yes" -a "x$AVX_CFLAGS" = x; then
		    AX_CHECK_COMPILER_FLAGS(-mavx, [AVX_CFLAGS="-mavx"],
			[AC_MSG_WARN([Need a version of gcc with -mavx])])
		fi

		if test "$have_altivec" = "yes" -a "x$ALTIVEC_CFLAGS" = x; then
		    # -DFAKE__VEC__ is a workaround because gcc-3.3 does not
		    # #define __VEC__ with -maltivec.
		    AX_CHECK_COMPILER_FLAGS(-faltivec, [ALTIVEC_CFLAGS="-faltivec"],
			[AX_CHECK_COMPILER_FLAGS(-maltivec -mabi=altivec,
			    [ALTIVEC_CFLAGS="-maltivec -mabi=altivec -DFAKE__VEC__"],
			    [AX_CHECK_COMPILER_FLAGS(-fvec, [ALTIVEC_CFLAGS="-fvec"],
				[AC_MSG_WARN([Need a version of gcc with -maltivec])])])])
		fi

		if test "$have_neon" = "yes" -a "x$NEON_CFLAGS" = x; then
		    AX_CHECK_COMPILER_FLAGS(-mfpu=neon, [NEON_CFLAGS="-mfpu=neon"],
			[AC_MSG_WARN([Need a version of gcc with -mfpu=neon])])
		fi
	esac
	AC_SUBST(SSE2_CFLAGS)
	AC_SUBST(PIPE_CFLAGS)
	AC_SUBST(TREE_VECTORIZE_CFLAGS)
	AC_SUBST(AVX_CFLAGS)
	AC_SUBST(ALTIVEC_CFLAGS)
	AC_SUBST(NEON_CFLAGS)
	OPTIMISATION_DETECTED="-O2 ${SSE2_CFLAGS} ${TREE_VECTORIZE_CFLAGS} ${AVX_CFLAGS} ${ALTIVEC_CFLAGS} ${NEON_CFLAGS} ${PIPE_CFLAGS}"
	AC_SUBST(OPTIMISATION_DETECTED)
else
	if test "$debug" = "yes"; then
		DEBUG_FLAGS+=" -g"
	fi
	if test "$profiling" = "yes"; then
		DEBUG_FLAGS+=" -pg"
	fi
fi
#VERBOSE_FLAGS="-DVERBOSE2"
#VERBOSE_FLAGS="-DVERBOSE3"
VERBOSE_FLAGS=""
if test "$verbose" = "yes"; then
	VERBOSE_FLAGS+=" -DVERBOSE"
fi
WERROR_FLAGS=""
if test "$werror" == "yes"; then
	WERROR_FLAGS+=" -Werror"
fi
if test "$valgrind" == "yes"; then
	OPTIMISATION_DETECTED=""
	if test "$debug" = "no"; then DEBUG_FLAGS+=" -g"; fi
fi

LDFLAGS="${MY_LDFLAGS} ${LDFLAGS} ${OPTIMISATION_DETECTED} ${DEBUG_FLAGS} ${VERBOSE_FLAGS} ${WERROR_FLAGS} ${HAVE_PTHREADS} ${PTHREAD_LIBS}"
CFLAGS="${MY_CFLAGS} ${CFLAGS} ${OPTIMISATION_DETECTED} ${DEBUG_FLAGS} ${VERBOSE_FLAGS} ${WERROR_FLAGS} ${HAVE_PTHREADS} ${PTHREAD_CFLAGS}"
CPPFLAGS="${MY_CPPFLAGS} ${CPPFLAGS} ${OPTIMISATION_DETECTED} ${DEBUG_FLAGS} ${VERBOSE_FLAGS} ${WERROR_FLAGS} ${HAVE_PTHREADS} ${PTHREAD_CFLAGS}"

echo "OPTIMISATIONS * USED *: ${OPTIMISATION_DETECTED}"
echo "VERBOSE FLAGS : ${VERBOSE_FLAGS}"
echo "DEBUG FLAGS : ${DEBUG_FLAGS}"
echo "VALGRIND enabled : $valgrind"
echo "CFLAGS USED: $CFLAGS"
echo "CPPFLAGS USED: $CPPFLAGS"

AM_COND_IF([DATA_FILE_CONSISTENCY_CHECKS],
[
	CFLAGS="${CFLAGS} -DDATA_FILE_CONSISTENCY_CHECKS=1"
	CPPFLAGS="${CPPFLAGS} -DDATA_FILE_CONSISTENCY_CHECKS=1"
],[echo]
)

AC_SUBST(LDFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_SUBST(AUTHOR, "Andreas Hadjiprocopis (livantes\@soi.city.ac.uk)")
AC_SUBST(BUILTDATE, `date`)
AC_SUBST(TITLE, "A lexical analyser and parser for Neural Network Scripts")
AC_SUBST(DATESTARTED, "Tue Mar 25 1997")
AC_SUBST(SEARCHPATH, $PATH)

AC_CONFIG_FILES([src/C/xa_stats.pc])

AC_OUTPUT(
	Makefile
	src/Makefile
	src/C/Makefile
)
