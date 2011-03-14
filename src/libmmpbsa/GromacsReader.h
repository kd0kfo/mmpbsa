#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define USE_GROMACS 1
#ifdef USE_GROMACS
#ifndef GROMACSREADER_H
#define GROMACSREADER_H


#include <string>
#include <valarray>
#include <vector>
#include <iostream>
#include <sstream>

#include "globals.h"
#include "mmpbsa_exceptions.h"

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



namespace mmpbsa_io{

typedef struct {

	std::map<std::string,size_t> indices;
}gromacs_fileindex;


typedef struct {
	size_t atom;//index of first atom of the molecule in the whole systems array list in SanderParm.

	//difference between gromacs index value and sanderparm.
	//For example, for bond values, index into sanderparm array = gromacs index - f_bonds
	size_t f_bonds,f_g96bonds,f_angles,f_g96angles,f_pdihs/*,f_idihs*/,f_lj14;
	size_t residue;

	//Offset for first array entry into sanderparm when putting values into the array *when
	//loading the molecule*. Suffix 'h' means with hydrogen array in sanderparm; Suffix 'a'
	//is the array without hydrogen.
	size_t bndh,bnda,thetah,theta,phih,phia,lj14/*Lennard jones 1-4*/;
}gromacs_idx_offsets;

/**
 * Reads Gromacs Trajectory Files (.trr)
 *
 * Loads coordinates cooresponding to the specified frame (snap shot) from the provided
 * file. If the requested frame is beyond the total number of frames, crds is set to an
 * empty array.
 */
void load_gmx_trr(const std::string& filename,std::valarray<mmpbsa_t>& crds,size_t frame_number,const size_t* natom_limit = 0);

/**
 * Determines whether or not the desired f rame is listed in the trr trajectory file.
 */
bool gmx_trr_eof(const std::string& filename,size_t frame_number);

/*
 *
 */
size_t total_gmx_trr_frames(const std::string& filename);

//void gmxfile_index(std::iostream& ascii_file);

std::vector<size_t> allowed_gmx_energies();
size_t& get_gmxarray_offset(mmpbsa_io::gromacs_idx_offsets& offsets,const size_t& gmx_bond_type);
size_t& get_gmxfunct_offset(mmpbsa_io::gromacs_idx_offsets& offsets,const size_t& gmx_bond_type);

}//end namespace mmpbsa_io

/**
 * Initializes the gromacs offset structure.
 */
void init(mmpbsa_io::gromacs_idx_offsets& offset);

/**
 * Initializes the gromacs offset structure, using void init(mmpbsa_io::gromacs_idx_offsets& offset)
 */
void init(mmpbsa_io::gromacs_idx_offsets* offset);

#endif//GROMACSREADER_H
#endif//use gromacs



