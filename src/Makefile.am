bin_PROGRAMS = mmpbsa
mmpbsa_CPPFLAGS = $(XML_CPPFLAGS) -I$(MEAD_PATH)/include/ -I./libmmpbsa $(BOINC_CPPFLAGS)
mmpbsa_LDFLAGS =  -L$(MEAD_PATH)/lib $(XML_LIBS) -L./libmmpbsa -L./molsurf $(BOINC_LDFLAGS)
mmpbsa_LDADD = -lmmpbsa -lmolsurf -lmead $(BOINC_LIBS)
mmpbsa_DEPENDENCIES = libmmpbsa
mmpbsa_SOURCES = mmpbsa.cpp
