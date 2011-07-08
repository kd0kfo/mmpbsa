bin_PROGRAMS = xyzr2sas
lib_LIBRARIES = libmolsurf.a

libmolsurf_adir = $(libdir)/lib
libmolsurf_a_SOURCES = molsurf.c
libmolsurf_a_includedir = $(includedir)/molsurf
libmolsurf_a_include_HEADERS = molsurf.h

xyzr2sas_SOURCES = xyzr2sas.c molsurf.c
xyzr2sas_LDADD = -lm
