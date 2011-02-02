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

void get_gromacs_forcefield(const char* fn,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list);

#if 0 // probably want to get rid of these
mmpbsa::SanderParm* gmxtpr2parmtop(const std::string& filename);
mmpbsa::SanderParm* gmxtpr2parmtop(std::iostream& gmxtop);
mmpbsa::SanderParm* gmxtpr2parmtop(const char* fn);
#endif

}//end namespace mmpbsa_io


#endif // TRAJECTORYREADER_H

