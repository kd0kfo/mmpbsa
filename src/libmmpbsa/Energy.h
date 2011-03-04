#ifndef MMPBSA_ENERGY_H
#define MMPBSA_ENERGY_H

#include <iostream>
#include <vector>
#include <valarray>

#include "structs.h"
#include "mmpbsa_utils.h"
#include "mmpbsa_exceptions.h"

namespace mmpbsa
{

	/**
	 * Given a list of bond_t structures and positions, the total bond energy is calculated
	 */
	mmpbsa_t bond_energy_calc(const std::vector<bond_t>& bonds,const std::valarray<mmpbsa_t>& crds);

	/**
	 * Given a list of angle_t structures and positions, the total angle energy is calculated
	 */
	mmpbsa_t angle_energy_calc(const std::vector<angle_t>& angles,const std::valarray<mmpbsa_t>& crds);

	/**
	 * Given a list of dihedral_t structures and positions, the total dihedral energy is calculated
	 */
	mmpbsa_t dihedral_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::valarray<mmpbsa_t>& crds);

	/**
	 * Given a list of dihedral_t structures and positions, the total Van der Waals energy between the
	 * first and fourth atoms is calculated
	 */
	mmpbsa_t vdw14_energy_calc(const std::vector<dihedral_t>& dihedrals,const std::valarray<mmpbsa_t>& crds,const mmpbsa_t& inv_scnb);

	/**
	 * Given a list of dihedral_t structures and positions, the total electrostatic energy between the
	 * first and fourth atoms is calculated
	 */
	mmpbsa_t elstat14_energy_calc(const std::vector<dihedral_t>& dihedrals, const std::vector<atom_t>& atoms, const std::valarray<mmpbsa_t>& crds, const mmpbsa_t& inv_scee, const mmpbsa_t& dielc);

	/**
	 * Given a list of atom_t structures and positions, the total Van der Waals energy is calculated
	 */
	mmpbsa_t vdwaals_energy(const std::vector<atom_t>& atoms, const std::vector<lj_params_t>& lj_params,const std::valarray<mmpbsa_t>& crds);

	/**
	 * Given a list of atom_t structures and positions, the total electrostatic energy is calculated
	 */
	mmpbsa_t total_elstat_energy(const std::vector<mmpbsa::atom_t>& atoms, const std::valarray<mmpbsa_t>& crds, const mmpbsa_t& coulomb_const = 1);
}

#endif//MMPBSA_ENERGY_H

