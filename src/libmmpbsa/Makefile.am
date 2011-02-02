lib_LIBRARIES = libmmpbsa.a
libmmpbsa_adir=$(libdir)
libmmpbsa_a_CPPFLAGS = -Wall  $(XML_CPPFLAGS) -I$(MEAD_PATH)/include/ -I../ $(BOINC_CPPFLAGS)
libmmpbsa_a_SOURCES = EmpEnerFun.cpp EMap.cpp EnergyInfo.cpp SanderInterface.cpp MeadInterface.cpp SanderParm.cpp mmpbsa_exceptions.cpp mmpbsa_utils_templates.cpp mmpbsa_utils.cpp XMLParser.cpp XMLNode.cpp mmpbsa_io.cpp StringTokenizer.cpp MMPBSAState.cpp Energy.cpp structs.cpp

libmmpbsa_a_includedir = $(includedir)/libmmpbsa
libmmpbsa_a_include_HEADERS = EmpEnerFun.h EMap.h EnergyInfo.h SanderInterface.h MeadInterface.h SanderParm.h mmpbsa_exceptions.h mmpbsa_utils.h mmpbsa_io.h StringTokenizer.h XMLParser.h XMLNode.h MMPBSAState.h Energy.h structs.h

if BUILD_WITH_GZIP
libmmpbsa_a_SOURCES += Zipper.cpp
libmmpbsa_a_include_HEADERS += Zipper.h
endif

if BUILD_WITH_GROMACS
libmmpbsa_a_SOURCES += FormatConverter.cpp GromacsReader.cpp
libmmpbsa_a_include_HEADERS += FormatConverter.h GromacsReader.h
libmmpbsa_a_CPPFLAGS += -I $(GROMACS_PATH)/include/
endif

