#ifndef MMPBSA_STRUCTS_H
#define MMPBSA_STRUCTS_H

#include "mmpbsa_utils.h"

namespace mmpbsa{

typedef struct {
	size_t name_index;///<Index into name array
	int atomic_number;
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




}//end namespace mmpbsa

#endif//MMPBSA_STRUCTS_H

