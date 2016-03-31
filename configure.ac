#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
# 5 argument version only available with aclocal >= 2.64
AC_INIT([CAM], [1.1.0], [], [cam], [http://www.bioinf.uni-freiburg.de/])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_SRCDIR([src/ITS_CSP.h])
AC_CANONICAL_HOST

# backwards compatibility for autoconf >= 2.64
# PACKAGE_URL should be the fifth argument of AC_INIT
m4_ifndef([AC_PACKAGE_URL],
	[m4_define([AC_PACKAGE_URL], [http://www.bioinf.uni-freiburg.de/])
	AC_DEFINE(PACKAGE_URL, "AC_PACKAGE_URL", [Package URL])
	AC_SUBST(PACKAGE_URL, AC_PACKAGE_URL)])

lt_enable_auto_import=""
case "$host_os" in 
	cygwin* | mingw* | cegcc*)
		AM_LDFLAGS="-Wl,--enable-auto-import $AM_LDFLAGS"
esac

# store current compiler flags to avoid default setup via AC_PROG_CXX and *_CC
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB # for VF2-library
# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS


AM_INIT_AUTOMAKE

AC_LANG([C++])


###############################################################################
# ANSI CONFORM CODE CHECK SETUP
###############################################################################

AC_MSG_CHECKING([whether to build strict ANSI C++ code])
ansi=no
AC_ARG_ENABLE([ANSI],
	[AS_HELP_STRING([--enable-ANSI],
	    [enable check for strict ANSI C++ source code (def=no)])],
  	[ansi="$enableval"])
AC_MSG_RESULT([$ansi])
if test x"$ansi" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -ansi -pedantic"
fi

###############################################################################
# SETUP FOR COMPILATION WITH 2011 ISO C++ STANDARD 
###############################################################################

AC_MSG_CHECKING([whether to build with 2011 ISO C++ standard])
cpp11=no
AC_ARG_ENABLE([c++11],
	[AS_HELP_STRING([--enable-c++11],
	    [enable build with 2011 ISO C++ standard (Check if supported by your compiler! def=no)])],
  	[cpp11="$enableval"])
AC_MSG_RESULT([$cpp11])
if test x"$cpp11" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -std=c++11"
fi


###############################################################################
# DEBUG SUPPORT SETUP
###############################################################################

AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	    [enable debug data generation (def=no)])],
  	[debuger="$enableval"])
AC_MSG_RESULT([$debuger])
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG])
	AM_CXXFLAGS="$AM_CXXFLAGS -g -O0 -Wall -Werror -Wno-uninitialized -Wno-deprecated"
else
	AC_DEFINE([NDEBUG])
	AM_CXXFLAGS="$AM_CXXFLAGS -O3 -fno-strict-aliasing -Wno-uninitialized -Wno-deprecated"
fi



###############################################################################
# GPROF PROFILER BUILT SUPPORT 
###############################################################################

AC_MSG_CHECKING([whether to build with '-pg' for gprof profiling])
profiling=no
AC_ARG_ENABLE([profiling],
	[AS_HELP_STRING([--enable-profiling],
	    [enable the built of gprof profiling information (-pg) (def=no)])],
  	[profiling="$enableval"])
AC_MSG_RESULT([$profiling])
if test x"$profiling" = x"yes"; then
	AM_CXXFLAGS="$AM_CXXFLAGS -pg"
fi



###############################################################################
# STATIC LINKING SUPPORT
###############################################################################

AC_ARG_ENABLE([static-linking],
    [AS_HELP_STRING([--enable-static-linking],
	    [enables static linking of the executables (default=disabled)])],
    [enable_static_linking=$enableval],
    [enable_static_linking=no])

if test "$enable_static_linking" = "yes"; then
	AM_LDFLAGS="-static $AM_LDFLAGS"
else
	AS_CASE([$host_os],
			[cygwin*],
			[AC_MSG_NOTICE([NOTE: '-fPIC' g++ option not used for Windows Cygwin build])],
			[AM_CXXFLAGS="$AM_CXXFLAGS -fPIC"])
fi
AC_MSG_CHECKING([whether to link only with static libraries])
AC_MSG_RESULT([$enable_static_linking])





###############################################################################
# BOOST LIBRARY PATH SUPPORT
###############################################################################

# BOOST header path support, if not installed in usual directories
AC_ARG_WITH([boost-include],
	[AS_HELP_STRING([--with-boost-include=prefix],
	    [try this for non-standard boost-header directory (exclude 'boost' sub-directory from prefix)])],
	[BOOSTPATHSET=1],
	[BOOSTPATHSET=0])
if test  "$BOOSTPATHSET" = "1" ; then
	AM_CXXFLAGS="-I$with_boost_include $AM_CXXFLAGS"
fi


##########################################################################
# adding the ggl library (e.g. with static name 'libggl.a')
##########################################################################

# adding the lib to the files to link
LIBS="-lgglchem -lggl -lsgm_vf2 -lsgm -lvf2 $LIBS"

# introduce the optional configure parameter for a non-standard install prefix of ggl
AC_ARG_WITH([ggl],
    [AS_HELP_STRING([--with-ggl=prefix],
        [try this for a non-standard install prefix of the ggl library])],
    [gglPATHSET=1],
    [gglPATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $gglPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_ggl/include"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_ggl/lib"
fi

##########################################################################
# adding the gecode library (e.g. with static name 'libgecode.a')
##########################################################################

# adding the lib to the files to link
LIBS="-lgecodesearch -lgecodeint -lgecodekernel -lgecodesupport $LIBS"
#LIBS="-lgecodesearch -lgecodeminimodel -lgecodeset -lgecodeint -lgecodekernel -lgecodesupport -lgecodedriver $LIBS"

# introduce the optional configure parameter for a non-standard install prefix of gecode
AC_ARG_WITH([gecode],
    [AS_HELP_STRING([--with-gecode=prefix],
        [try this for a non-standard install prefix of the gecode library])],
    [gecodePATHSET=1],
    [gecodePATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $gecodePATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_gecode/include -pthread"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_gecode/lib -pthread"
fi

##########################################################################

###############################################################################
# DEPENDENCY CHECKS
###############################################################################


###########################################
# check for stl headers
AC_CHECK_HEADERS([vector])

###########################################
# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T



##########################################################################
# check for ggl headers
##########################################################################

# store current CXXFLAGS and merge with AM_CXXFLAGS for compilation check   
OLD_CXXFLAGS=$CXXFLAGS;
 # merge compiler flags and ignore warnings for checking purpose
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS -w"
# check for ggl library headers   
AC_MSG_CHECKING([for the Graph Grammar library headers])
# try to compile a file that includes a header of the library ggl
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[ #include <sgm/HashMap.hh> ]],[[]])],
    [AC_MSG_RESULT([found])
        FOUND_ggl=1;],
    [AC_MSG_RESULT([not found])
        FOUND_ggl=0;])

# reset original CXXFLAGS
CXXFLAGS=$OLD_CXXFLAGS

###########################################
# check for boost library headers   
OLD_CXXFLAGS=$CXXFLAGS;
 # merge compiler flags and ignore warnings for checking purpose
CXXFLAGS="$CXXFLAGS $AM_CXXFLAGS -w"
AC_MSG_CHECKING([for the Boost BGL C++ library headers])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[ #include <boost/graph/adjacency_list.hpp> ]],[[]])],
     [AC_MSG_RESULT([found])
		FOUND_BOOST=1;],
     [AC_MSG_RESULT([not found])
		FOUND_BOOST=0;])
CXXFLAGS=$OLD_CXXFLAGS


###############################################################################
# COMPILER AND LINKER FLAG DISTRIBUTION
###############################################################################

# distribute the configure based compiler and linker setups
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([LIBS])
AC_SUBST([PREFIX_PATH], [${prefix}])


###############################################################################
# DEPENDENCY CHECK RESULT HANDLING
###############################################################################

STOPCONFIGURE=0;

if test $FOUND_ggl != 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([The Graph Grammar Library (GGL) was not found!])
	if test  $gglPATHSET = 1 ; then
		AC_MSG_NOTICE([ The path '$with_ggl' was set. Maybe wrong ???])
	else
		AC_MSG_NOTICE([ No non-default path was set. You might want to use '--with-ggl=PREFIX' ?!?])
	fi
	AC_MSG_NOTICE([])
	STOPCONFIGURE=1;
fi

if test $FOUND_BOOST != 1; then
	AC_MSG_NOTICE([])
	AC_MSG_NOTICE([The Boost Graph Library (BGL) was not found!])
	if test  $BOOSTPATHSET = 1 ; then
		AC_MSG_NOTICE([ The path '$with_boost_include' was set. Maybe wrong ???])
	else
		AC_MSG_NOTICE([ No non-default path was set. You might want to use '--with-boost-include=PREFIX' ?!?])
	fi
	AC_MSG_NOTICE([])
	STOPCONFIGURE=1;
fi

if test $STOPCONFIGURE != 0; then
	AC_MSG_ERROR([Essential dependencies are missing : cannot compile and stop here !])
fi


###############################################################################
# FINAL SETUP
###############################################################################


AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_OUTPUT


###############################################################################
# SETUP OUTPUT
###############################################################################

echo ""
echo "----------------------------------"
echo " Build Settings for ${PACKAGE_TARNAME} v${PACKAGE_VERSION}"
echo "----------------------------------"
echo "  host type            = $host_os"
echo "  install dir          = $prefix"
echo "  CC                   = $CC"
echo "  CFLAGS               = $AM_CPPFLAGS $CFLAGS"
echo "  CXX                  = $CXX"
echo "  CXXFLAGS             = $AM_CXXFLAGS $CXXFLAGS"
echo "  LDFLAGS              = $AM_LDFLAGS $LDFLAGS"
echo "  LIBS                 = $LIBS"
echo "----------------------------------"
echo ""
echo " You can run ./make now."
echo ""
