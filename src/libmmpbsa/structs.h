#ifndef MMPBSA_STRUCTS_H
#define MMPBSA_STRUCTS_H

#include <iostream>
#include <vector>
#include <string>
#include <set>

#include "globals.h"

namespace mmpbsa{

/**
 * Lennard Jones Coefficients, using kcal and Angstrom energy units.
 */
typedef struct{
	mmpbsa_t c6,c12;
}lj_params_t;

/**
 * Atom data structure.
 */
typedef struct {
	std::string name;///<Index into name array
	int atomic_number;
	mmpbsa_t charge;
	size_t atom_type;///<For Lennard Jones Parameters, etc
	std::set<size_t> exclusion_list;
}atom_t;

/**
 * Bond data structure
 */
typedef struct {
	mmpbsa_t energy_const;///<Proportionality constant for bond energy, in kcal.
	mmpbsa_t eq_distance;///<Equilibrium distance, in Angstroms.
}bond_energy_t;

/**
 * Structure of a bond pair.
 */
typedef struct {
	size_t atom_i,atom_j;///<Indices of atoms for bond pairs
	bond_energy_t* bond_energy;///<Pointer to bond's energy data.
}bond_t;

/**
 * Structure for an angle triplet
 */
typedef struct {
	size_t atom_i,atom_j,atom_k;///<Indices of atoms of angle triplets
	bond_energy_t* angle_energy;///<Pointer to bond's energy data.
}angle_t;

/**
 * Dihedral energy struture.
 */
typedef struct{
	mmpbsa_t periodicity;
	mmpbsa_t energy_const;///<Dihedral energy in kcal
	mmpbsa_t phase;///<radians
}dihedral_energy_t;

/**
 * Structure for dihedral triplet
 */
typedef struct {
	size_t atom_i,atom_j,atom_k,atom_l;///<Indices of atoms of angle triplets
	struct{mmpbsa_t c6,c12;}lj;///<Lennard Jones Parameters
	struct{bool is_improper,should_ignore_end_grp;}nonbonded_masks;
	dihedral_energy_t* dihedral_energy;///<Pointer to bond's energy data.
}dihedral_t;

/**
 * Forcefield structure.
 */
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

namespace mmpbsa_io{
/**
 * Trajectory data structure. Abstracts away the type of trajectory
 * file.
 *
 * This structure maintains the information needs for moving through
 * trajectory fields and retrieving data. When using the data
 * structure, the user does not need to know what type of trajectory
 * is being used. Methods operating on the structure handle the
 * file type accordingly.
 */
typedef struct {
	size_t curr_snap;
	std::iostream::streampos curr_pos;
	//for sander trajectories
	std::string* sander_filename;
	std::iostream* sander_crd_stream;
	size_t natoms;
	int ifbox;

	//for gromacs trajectories
	std::string* gromacs_filename;
	size_t num_gmx_frames;
}trajectory_t;

}//mmpbsa_io namespace

/**
 * Initializes the forcefield.
 */
void init(mmpbsa::forcefield_t* ff);

/**
 * Destroys the forcefield.
 *
 * This should be called for every force field to prevent memory leaks.
 */
void destroy(mmpbsa::forcefield_t* ff);


#endif//MMPBSA_STRUCTS_H

