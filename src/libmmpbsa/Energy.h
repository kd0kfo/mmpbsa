#ifndef MMPBSA_ENERGY_H
#define MMPBSA_ENERGY_H

#include <iostream>

#include "structs.h"

namespace mmpbsa
{

	mmpbsa_t bond_energy_calc(const std::vector<bond_t>& bonds,const std::valarray<mmpbsa_t>& crds);

	mmpbsa_t angle_energy_calc(const std::vector<angle_t>& angles,const std::valarray<mmpbsa_t>& crds);

	mmpbsa_t dihedral_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::valarray<mmpbsa_t>& crds);

	mmpbsa_t vdw14_energy_calc(const std::vector<dihedral_t>& dihedrals,const std::valarray<mmpbsa_t>& crds,const mmpbsa_t& inv_scnb);

	mmpbsa_t elstat14_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::vector<atom_t>& atoms, const std::valarray<mmpbsa_t>& crds, const mmpbsa_t& inv_scee, const mmpbsa_t& dielc);

	mmpbsa_t vdwaals_energy(const std::vector<atom_t>& atoms, const std::vector<lj_params_t>& lj_params,const std::valarray<mmpbsa_t>& crds);
}

#endif//MMPBSA_ENERGY_H

