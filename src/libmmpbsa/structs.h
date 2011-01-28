#ifndef MMPBSA_STRUCTS_H
#define MMPBSA_STRUCTS_H

#include <string>
#include <set>

#include "mmpbsa_utils.h"

namespace mmpbsa{

typedef struct{
	mmpbsa_t c6,c12;
}lj_params_t;

typedef struct {
	std::string name;///<Index into name array
	int atomic_number;
	mmpbsa_t charge;
	size_t atom_type;///<For Lennard Jones Parameters, etc
	std::set<size_t> exclusion_list;
}atom_t;

typedef struct {
	mmpbsa_t energy_const;///<Proportionality constant for bond energy.
	mmpbsa_t eq_distance;///<Equilibrium distance.
}bond_energy_t;

typedef struct {
	size_t atom_i,atom_j;///<Indices of atoms for bond pairs
	bond_energy_t* bond_energy;///<Pointer to bond's energy data.
}bond_t;

typedef struct {
	size_t atom_i,atom_j,atom_k;///<Indices of atoms of angle triplets
	bond_energy_t* angle_energy;///<Pointer to bond's energy data.
}angle_t;

typedef struct{
	mmpbsa_t periodicity;
	mmpbsa_t energy_const;
	mmpbsa_t phase;///<radians
}dihedral_energy_t;

typedef struct {
	size_t atom_i,atom_j,atom_k,atom_l;///<Indices of atoms of angle triplets
	struct{mmpbsa_t c6,c12;}lj;///<Lennard Jones Parameters
	struct{bool is_improper,should_ignore_end_grp;}nonbonded_masks;
	dihedral_energy_t* dihedral_energy;///<Pointer to bond's energy data.
}dihedral_t;


typedef struct {
	//indices
	std::vector<bond_t> bonds_with_H,bonds_without_H;
	std::vector<angle_t> angles_with_H,angles_without_H;
	std::vector<dihedral_t> dihedrals_with_H,dihedrals_without_H;

	//Data
	mmpbsa::bond_energy_t *bond_energy_data,*angle_energy_data;
	mmpbsa::dihedral_energy_t *dihedral_energy_data;
	std::vector<mmpbsa::lj_params_t> lj_params;

	//Constants
	mmpbsa_t inv_scnb, inv_scee, dielc, coulomb_const;
}forcefield_t;

}//end namespace mmpbsa

void init(mmpbsa::forcefield_t* ff);
void destroy(mmpbsa::forcefield_t* ff);


#endif//MMPBSA_STRUCTS_H

