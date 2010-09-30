#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT([mmpbsa], [0.6.7], [David.Coss@stjude.org])
AC_CONFIG_SRCDIR([src/libmmpbsa/])
AC_CONFIG_SRCDIR([src/])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

AH_TEMPLATE([USE_BOINC],[Define to use the BOINC library])
AH_TEMPLATE([USE_GRAPHICS],[Define to create the graphics application. Requires BOINC and GL])
AH_TEMPLATE([USE_JPEGLIB],[Define to use the JPEG library])
AH_TEMPLATE([WITHOUT_MOLSURF],[Define to prevent the use of molsurf, i.e. on Mac without proper signal.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

AC_SUBST(CUSTOM_CPPFLAGS,"-Wall")

AC_ARG_ENABLE([static],AS_HELP_STRING([--enable-static],[Compile Staticahlly]),
[AC_MSG_RESULT(Compiling MMPBSA statically)
AC_SUBST(CUSTOM_LDFLAGS," -static ")
AC_SUBST(CUSTOM_LIBS,"-lpthread")],)

AC_ARG_ENABLE([molsurf],AS_HELP_STRING([--enable-molsurf],[Compile and link to molsurf. Without molsurf, Surface Area calculations are not done.]),molsurf=$enableval)
if test x"$molsurf" == xno; then
	AC_DEFINE([WITHOUT_MOLSURF])
fi

AC_ARG_WITH([mead],AS_HELP_STRING([--with-mead@<:@=dir@:>@],[Include MEAD library @<:@Default=PREFIX/include@:>@]),with_mead=${withval},AC_MSG_RESULT(Assuming mead is installed in system paths))
if test x$with_mead != x; then
   AC_SUBST(MEAD_PATH,$with_mead)
fi

AC_ARG_WITH([boinc],AS_HELP_STRING([--with-boinc@<:@=dir@:>@],[Include boinc under specified path @<:@Default=no@:>@]),with_boinc=${withval},AC_MSG_RESULT(Not using BOINC))
if test x$with_boinc != x; then
   AC_SUBST(BOINC_CPPFLAGS,"-I$with_boinc/include/")
   AC_SUBST(BOINC_LIBS,"-lboinc_api -lboinc")
   AC_SUBST(BOINC_LDFLAGS,"-L$with_boinc/lib")
   AC_DEFINE([USE_BOINC],[1])   
fi

AC_ARG_WITH([graphics],AS_HELP_STRING([--with-graphics@<:@=dir@:>@],[Include the graphics app. Requires BOINC and GL. @<:@Default=no@:>@]),
[AC_MSG_RESULT(Using Graphics $with_graphics)
AC_SUBST(GRAPHICS_CPPFLAGS,"-DAPP_GRAPHICS -I$with_graphics")
AC_SUBST(GRAPHICS_LIBS,"-lboinc_graphics2 -lopengl32 -lglu32 -lgdi32")
AC_SUBST(GRAPHICS_LDFLAGS,"-L/usr/X11R6/lib")
AC_DEFINE([USE_GRAPHICS],[1])
 ],
AC_MSG_RESULT(Not creating a graphics application))

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([floor pow sqrt strlcpy])

AC_CONFIG_FILES([Makefile
		 src/libmmpbsa/Makefile
                 src/Makefile
                 src/molsurf/Makefile])
AC_OUTPUT
