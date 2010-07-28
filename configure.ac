#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT([mmpbsa], [0.2], [David.Coss@stjude.org])
AC_CONFIG_SRCDIR([src/libmmpbsa/])
AC_CONFIG_SRCDIR([src/])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

AH_TEMPLATE([__USE_BOINC__],[Define to use the BOINC library])
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

AM_PATH_XML2(2.6.0)

AC_ARG_WITH([mead],AS_HELP_STRING([--with-mead],[Include MEAD library]),with_mead=${withval},AC_MSG_RESULT(Assuming mead is installed in system paths))
if test x$with_mead != x; then
   AC_SUBST(MEAD_PATH,$with_mead)
fi

AC_ARG_WITH([boinc],AS_HELP_STRING([--with-boinc],[Include boinc under specified path]),with_boinc=${withval},AC_MSG_RESULT(Not using BOINC))
if test x$with_boinc != x; then
   AC_SUBST(BOINC_CPPFLAGS,"-I$with_boinc/client -I$with_boinc/api -I$with_boinc/lib -I$with_boinc")
   AC_SUBST(BOINC_LIBS,"-lboinc_api -lboinc -lpthread")
   AC_SUBST(BOINC_LDFLAGS,"-L$with_boinc/api -L$with_boinc/lib")
   AC_DEFINE([__USE_BOINC__],[1])   
fi

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile
		 src/libmmpbsa/Makefile
                 src/Makefile
                 src/molsurf/Makefile])
AC_OUTPUT
