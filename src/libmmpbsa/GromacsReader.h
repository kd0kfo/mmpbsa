#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_GROMACS
#ifndef GROMACSREADER_H

#include "mmpbsa_utils.h"

#include <string>
#include <valarray>
#include <iostream>

#ifndef eCPP_OK
//gromacs stuff
#ifdef HAVE_CONFIG_H
#include "gromacs/config.h"
#endif

//#include "gromacs/macros.h"
#include "gromacs/futil.h"
#include "gromacs/statutil.h"
#include "gromacs/sysstuff.h"
#include "gromacs/txtdump.h"
#include "gromacs/gmx_fatal.h"
#include "gromacs/xtcio.h"
#include "gromacs/enxio.h"
#include "gromacs/smalloc.h"
#include "gromacs/names.h"
#include "gromacs/gmxfio.h"
#include "gromacs/tpxio.h"
#include "gromacs/trnio.h"
#include "gromacs/txtdump.h"
//#include "gromacs/gmxcpp.h"
#include "gromacs/checkpoint.h"
#include "gromacs/mtop_util.h"
#include "gromacs/sparsematrix.h"
#include "gromacs/mtxio.h"
#endif//gromacs stuff

#if 0
#include <math.h>
#include "main.h"
#endif


namespace mmpbsa_io{

typedef struct {

	std::map<std::string,size_t> indices;
}gromacs_fileindex;


/**
 * Reads Gromacs Trajectory Files (.trr)
 *
 * Loads coordinates cooresponding to the specified frame (snap shot) from the provided
 * file. If the requested frame is beyond the total number of frames, crds is set to an
 * empty array.
 */
void load_gmx_trr(const std::string& filename,std::valarray<mmpbsa_t>& crds,size_t frame_number);

 gmxfile_index(std::iostream& ascii_file);

}//end namespace mmpbsa_io

#endif//GROMACSREADER_H
#endif//use gromacs



