bin_PROGRAMS = mmpbsa mmpbsa_graphics
mmpbsa_CPPFLAGS = $(XML_CPPFLAGS) -I$(MEAD_PATH)/include/ -I./libmmpbsa $(BOINC_CPPFLAGS) $(GRAPHICS_CPPFLAGS)
mmpbsa_LDFLAGS =  -L$(MEAD_PATH)/lib $(XML_LIBS) -L./libmmpbsa -L./molsurf $(BOINC_LDFLAGS) $(GRAPHICS_LDFLAGS)
mmpbsa_LDADD = -lmmpbsa -lmolsurf -lmead $(BOINC_LIBS) $(GRAPHICS_FLAGS)
mmpbsa_DEPENDENCIES = libmmpbsa
mmpbsa_SOURCES = mmpbsa.cpp

mmpbsa_graphics_CPPFLAGS = $(BOINC_CPPFLAGS) $(GRAPHICS_CPPFLAGS)
mmpbsa_graphics_LDFLAGS = $(BOINC_LDFLAGS) $(GRAPHICS_LDFLAGS)
mmpbsa_graphics_LDADD = $(BOINC_LIBS) $(GRAPHICS_LIBS)
mmpbsa_graphics_SOURCES = mmpbsa_graphics.cpp
