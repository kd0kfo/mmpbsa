#ifndef FORMATCONVERTER_H
#define FORMATCONVERTER_H

#include "globals.h"
#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"
#include "structs.h"
#include "GromacsReader.h"
#include "MMPBSAState.h"

//gromacs stuff
#ifdef USE_GROMACS
#include "gromacs/types/idef.h"
#include "gromacs/types/atoms.h"
#endif

namespace mmpbsa_io{

/**
 * Sets up the forcefield data from a Gromacs (http://gromacs.org) .tpr
 * topology file, which has a filename stored in fn.
 *
 * split_ff should be a pointer to an array mmpbsa::forcefield_t[3], which
 * will then store the complex, receptor and ligand in the array, in that order.
 * atom_lists should be a pointer of vectors which will be filled with atom
 * data in a manner similar to split_ff.
 *
 * mol_list will be filled with mmpbsa::MMPBSAState::MOLECULE values that indicate
 * whether atoms belong to the receptor or ligand. This is to be used later to
 * match coordinates, which are stored a vector, to their atoms.
 *
 * The optional paraters for receptor and ligand positions give the users the
 * option to specify which molecules are the receptor and ligand. This should be a
 * zero-indexed list corresponding to the molecule list in the .tpr file.
 * If no lists are provide, which is the default case, it is assumed that
 * the first molecule is the receptor and all other molecules not named
 * "SOL" are ligand.
 */
void get_gromacs_forcefield(const char* fn,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list,const std::set<size_t>* receptor_pos = 0, const std::set<size_t>* ligand_pos = 0);

}//end namespace mmpbsa_io

/**
 * Outputs the forcefield data to the provided file.
 */
void dump_forcefield(std::vector<mmpbsa::atom_t>* atoms, mmpbsa::forcefield_t* ff,const char* filename);

#endif // TRAJECTORYREADER_H

