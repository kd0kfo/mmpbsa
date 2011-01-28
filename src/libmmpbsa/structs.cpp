#include "structs.h"


void init(mmpbsa::forcefield_t* ff)
{
	ff->angle_energy_data = 0;
	ff->bond_energy_data = 0;
	ff->dihedral_energy_data = 0;
	ff->coulomb_const = 1;
	ff->inv_scee = 1;
	ff->inv_scnb = 1;
	ff->dielc = 1;
#if 0//not default here
	ff->inv_scee = 1/DEFAULT_SCNB;
	ff->inv_scnb = 1/DEFAULT_SCEE;
	ff->dielc = DEFAULT_DIELC;
#endif
}

void destroy(mmpbsa::forcefield_t* ff)
{
	delete [] ff->angle_energy_data;
	delete [] ff->bond_energy_data;
	delete [] ff->dihedral_energy_data;

	ff->angle_energy_data = 0;
	ff->bond_energy_data = 0;
	ff->dihedral_energy_data = 0;
}
