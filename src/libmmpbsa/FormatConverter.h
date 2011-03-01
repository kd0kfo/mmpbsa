#ifndef FORMATCONVERTER_H
#define FORMATCONVERTER_H

#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"
#include "mmpbsa_utils.h"
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
 * Sets up the forcefield data.
 *
 * mol_list is created as a mask for atom_lists, as to whether an atom is in the Receptor or ligand
 *
 * If the receptor and ligand starting position sets are null or empty,
 * it is assumed that the first molblock is the receptor and the reset are ligand, except those named SOL.
 *
 */
void get_gromacs_forcefield(const char* fn,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list,const std::set<size_t>* receptor_pos = 0, const std::set<size_t>* ligand_pos = 0);

}//end namespace mmpbsa_io

void dump_forcefield(std::vector<mmpbsa::atom_t>* atoms, mmpbsa::forcefield_t* ff,const char* filename);

#endif // TRAJECTORYREADER_H

