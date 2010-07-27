bin_PROGRAMS = mmpbsa
mmpbsa_CPPFLAGS = $(XML_CPPFLAGS) -I$(MEAD_PATH)/include/ -I./libmmpbsa
mmpbsa_LDFLAGS =  -L$(MEAD_PATH)/lib $(XML_LIBS) -L./libmmpbsa -L./molsurf
mmpbsa_LDADD = -lmmpbsa -lmolsurf -lmead
mmpbsa_DEPENDENCIES = libmmpbsa
mmpbsa_SOURCES = mmpbsa.cpp
