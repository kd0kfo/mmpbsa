#include "FormatConverter.h"

bool mmpbsa_io::get_next_snap(mmpbsa_io::trajectory_t& traj, std::valarray<mmpbsa_t>& snapshot)
{

#ifdef USE_GROMACS
	if(traj.gromacs_filename != 0)
	{
		mmpbsa_io::load_gmx_trr(*traj.gromacs_filename,snapshot,traj.curr_snap++);
		return snapshot.size() != 0;
	}
#endif

	if(traj.natoms == 0 || traj.sander_crd_stream == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::get_next_snap: Sander parameters and/or sander coordinate stream is missing.",mmpbsa::DATA_FORMAT_ERROR);
	return get_next_snap(*traj.sander_crd_stream,snapshot,traj.natoms,(traj.ifbox > 0));
}


void mmpbsa_io::seek(mmpbsa_io::trajectory_t& traj,const size_t& snap_pos)
{
	size_t i = traj.curr_snap;
	if(traj.sander_crd_stream != 0)
	{
		if(traj.natoms == 0)
			throw mmpbsa::MMPBSAException("mmpbsa_io::seek: Trajectory cannot be read without paramters. However, sander paramter object is a null pointer.",mmpbsa::NULL_POINTER);
		bool isPeriodic = (traj.ifbox > 0);//Are periodic boundary conditions used?
		try
		{
			for(;i<snap_pos;i++)
			{
				mmpbsa_io::skip_next_snap(*traj.sander_crd_stream,traj.natoms,isPeriodic);

			}//after this for loop, the trajFile is pointing to the beginning of currState.currentSnap
		}
		catch(mmpbsa::MMPBSAException e)
		{
			std::ostringstream error;
			error << "mmpbsa_io::seek: Died reading snap shots on snap number " << i
					<< std::endl << " Error message: " <<  e.what();
			throw mmpbsa::MMPBSAException(error,e.getErrType());
		}
	}
	traj.curr_snap = snap_pos;
}




void mmpbsa_io::default_trajectory(mmpbsa_io::trajectory_t& traj)
{
	traj.sander_crd_stream = 0;
	traj.natoms = 0;
	traj.ifbox = 0;

	traj.gromacs_filename = 0;
	traj.curr_snap = 0;
}

void mmpbsa_io::destroy_trajectory(mmpbsa_io::trajectory_t& traj)
{
	delete traj.sander_crd_stream;
	delete traj.gromacs_filename;
}

mmpbsa_io::trajectory_t mmpbsa_io::open_trajectory(const std::string& filename)
{
	trajectory_t returnMe;
	bool is_sander = true;
	mmpbsa_io::default_trajectory(returnMe);

#ifdef USE_GROMACS
	if(filename.find(".trr") != std::string::npos)
	{
		returnMe.gromacs_filename = new std::string(filename);
		return returnMe;
	}
#endif
	std::fstream trajDiskFile(filename.c_str(),std::ios::in);
	if(!trajDiskFile.good())
		throw mmpbsa::MMPBSAException("mmpbsa_io::open_trajectory: Unable to read from trajectory file",mmpbsa::BROKEN_TRAJECTORY_FILE);
	returnMe.sander_crd_stream = new std::stringstream;
	mmpbsa_io::smart_read(*returnMe.sander_crd_stream,trajDiskFile,&filename);
	trajDiskFile.close();
	return returnMe;

}

bool mmpbsa_io::eof(trajectory_t& traj)
{
	if(traj.sander_crd_stream == 0)
	{
#ifdef USE_GROMACS
		if(traj.gromacs_filename != 0)
			return mmpbsa_io::gmx_trr_eof(*traj.gromacs_filename,traj.curr_snap);
#endif
		return false;
	}
	return traj.sander_crd_stream->eof();
}


std::string mmpbsa_io::get_traj_title(mmpbsa_io::trajectory_t& traj)
{
	if(traj.sander_crd_stream == 0 && traj.gromacs_filename == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::get_traj_title: no trajectory provided.",mmpbsa::DATA_FORMAT_ERROR);
	if(traj.sander_crd_stream != 0)
		return mmpbsa_io::get_traj_title(*traj.sander_crd_stream);
	return *traj.gromacs_filename;
}

size_t gromacs_num_bond_types(const std::map<t_functype,size_t>& function_sizes)
{
	std::map<t_functype,size_t>::const_iterator function;
	size_t returnMe = 0;
	for(function = function_sizes.begin();function != function_sizes.end();function++)
	{
		switch(function->first)
		{
		case F_BONDS:
			returnMe += function->second;
		default:
			continue;
		}
	}
	return returnMe;
}

size_t gromacs_num_angle_types(const std::map<t_functype,size_t>& function_sizes)
{
	std::map<t_functype,size_t>::const_iterator function;
	size_t returnMe = 0;
	for(function = function_sizes.begin();function != function_sizes.end();function++)
	{
		switch(function->first)
		{
		case F_ANGLES:
			returnMe += function->second;
		default:
			continue;
		}
	}
	return returnMe;
}

size_t gromacs_num_dihedral_types(const std::map<t_functype,size_t>& function_sizes)
{
	std::map<t_functype,size_t>::const_iterator function;
	size_t returnMe = 0;
	for(function = function_sizes.begin();function != function_sizes.end();function++)
	{
		switch(function->first)
		{
		case F_PDIHS:case F_PIDIHS:
			returnMe += function->second;
		default:
			continue;
		}
	}
	return returnMe;
}

size_t gromacs_num_lj_parameters(const std::map<t_functype,size_t>& function_sizes)
{
	std::map<t_functype,size_t>::const_iterator function;
	size_t returnMe = 0;
	for(function = function_sizes.begin();function != function_sizes.end();function++)
	{
		switch(function->first)
		{
		case F_LJ:
			returnMe += function->second;
		default:
			continue;
		}
	}
	return returnMe;
}

void gromacs_adjust_start_indices(std::map<t_functype,size_t>& start_indices,const std::map<t_functype,size_t>& function_sizes)
{
	std::map<t_functype,size_t>::iterator index = start_indices.begin();
	for(;index != start_indices.end();index++)
	{
		switch(index->first)
		{
		case F_PIDIHS:
		{
			if(function_sizes.find(F_PDIHS) != function_sizes.end())
				index->second -= function_sizes.find(F_PDIHS)->second;//F_PDIHS before F_PIDIHS
			break;
		}
		}
	}
}

/**
 * Sets up the forcefield data.
 *
 * mol_list is created as a mask for atom_lists, as to whether an atom is in the Receptor or ligand
 *
 * It is assumed that the first molblock is the receptor and the rest are ligand, except for those named "SOL"
 * This should be changed.
 */
void mmpbsa_io::get_gromacs_forcefield(const char* fn,mmpbsa::forcefield_t** split_ff,std::vector<mmpbsa::atom_t>** atom_lists, std::valarray<mmpbsa::MMPBSAState::MOLECULE>& mol_list)
{
	using namespace mmpbsa;
	mmpbsa_io::gromacs_idx_offsets offsets;
	size_t atom_offset = 0;
	//unused
//	FILE *gp;
	int         fp,atot;
//	t_state     state;
	rvec        *f=NULL;
//	t_inputrec  ir;
	t_tpxheader tpx;
//	t_topology top;
	gmx_mtop_t  mtop;
//	gmx_groups_t *groups;
//	gmx_ffparams_t* force_field;
	const char* mdpfn = 0;
	bool bSysTop,bShowNumbers;
	bSysTop = false;
	bShowNumbers = true;

	static const mmpbsa_t joules2cal = 0.23889;
	static const mmpbsa_t angst2nm_sqrd = 0.01;
	static const mmpbsa_t nm2angst_6 = 1e+6;
	static const mmpbsa_t nm2angst_12 = 1e+12;
 	static const mmpbsa_t nm2angst = 10;
 	static const mmpbsa_t charge_units = 18.2182634799;

	if(split_ff == 0 || atom_lists == 0)
		throw mmpbsa::MMPBSAException("mmpbsa::MMPBSAException: Null pointer supplied for atom list and/or force field.",mmpbsa::NULL_POINTER);

	*split_ff = new mmpbsa::forcefield_t[MMPBSAState::END_OF_MOLECULES];
	*atom_lists = new std::vector<mmpbsa::atom_t>[MMPBSAState::END_OF_MOLECULES];

	//for simplicity alias these here.
	forcefield_t& complex = split_ff[0][MMPBSAState::COMPLEX];
	forcefield_t& receptor = split_ff[0][MMPBSAState::RECEPTOR];
	forcefield_t& ligand = split_ff[0][MMPBSAState::LIGAND];
	for(size_t i = 0;i<MMPBSAState::END_OF_MOLECULES;i++)
	{
		init(&split_ff[0][i]);
	}


	if(fn == 0 || strlen(fn) == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::gmxtpr2parmtop: Filename is required.",mmpbsa::DATA_FORMAT_ERROR);

	printf("Make sure to check DIELECTRIC CONSTANTS AND GROMACS UNITS!!!!");
	//const std::vector<size_t>& interaction_types = mmpbsa_io::allowed_gmx_energies();
	std::vector<size_t>::const_iterator interaction_type;

	std::map<t_functype,size_t> function_start_index,function_sizes;

	read_tpx(fn,NULL,NULL,&atot,0,0,0,&mtop);
	if (!mdpfn) {
		//top = gmx_mtop_t_to_t_topology(&mtop);
		if (available(stdout,&tpx,0,fn)) {
			//setup constants
			init(offsets);

			//Determine beginning of function type data. This will be used to
			//create dynamic arrays of interaction data.
			for(size_t i = 0;i<mtop.ffparams.ntypes;i++)//get offset values
			{
				if(function_sizes.find(mtop.ffparams.functype[i]) == function_sizes.end())
				{
					function_sizes[mtop.ffparams.functype[i]] = 1;
					function_start_index[mtop.ffparams.functype[i]] = i;
				}
				else
					function_sizes[mtop.ffparams.functype[i]]++;

			}

			//Allocate memory of energy data.
			complex.bond_energy_data = new bond_energy_t[gromacs_num_bond_types(function_sizes)];
			complex.angle_energy_data = new bond_energy_t[gromacs_num_angle_types(function_sizes)];
			complex.dihedral_energy_data = new dihedral_energy_t[gromacs_num_dihedral_types(function_sizes)];
			complex.lj_params.reserve(gromacs_num_lj_parameters(function_sizes));

			//Adjust start indices so that bond types may be placed in above allocated arrays
			//with similar bond types without overlapping.
			gromacs_adjust_start_indices(function_start_index,function_sizes);

			//Load energy data
			for(size_t i = 0;i<mtop.ffparams.ntypes;i++)//get offset values
			{
				switch(mtop.ffparams.functype[i])
				{
				case F_BONDS:
				{
					bond_energy_t new_energy;
					new_energy.energy_const = 0.5*mtop.ffparams.iparams[i].harmonic.krA*joules2cal*angst2nm_sqrd;
					new_energy.eq_distance = mtop.ffparams.iparams[i].harmonic.rA*nm2angst;
					complex.bond_energy_data[i-function_start_index[mtop.ffparams.functype[i]]] = new_energy;
					break;
				}
				case F_ANGLES:
				{	bond_energy_t new_energy;
					new_energy.energy_const = 0.5*mtop.ffparams.iparams[i].harmonic.krA*joules2cal;
					new_energy.eq_distance = MMPBSA_DEG_TO_RAD*mtop.ffparams.iparams[i].harmonic.rA;
					complex.angle_energy_data[i-function_start_index[mtop.ffparams.functype[i]]] = new_energy;
					break;
				}
				case F_PDIHS: case F_PIDIHS:
				{
					dihedral_energy_t new_energy;
					new_energy.energy_const = joules2cal*mtop.ffparams.iparams[i].pdihs.cpA;
					new_energy.periodicity = mtop.ffparams.iparams[i].pdihs.mult;
					new_energy.phase = MMPBSA_DEG_TO_RAD*mtop.ffparams.iparams[i].pdihs.phiA;
					complex.dihedral_energy_data[i-function_start_index[mtop.ffparams.functype[i]]] = new_energy;
					break;
				}
				case F_LJ:
					lj_params_t lj_params;
					lj_params.c12 = nm2angst_12*joules2cal*mtop.ffparams.iparams[i].lj.c12;
					lj_params.c6 = nm2angst_6*joules2cal*mtop.ffparams.iparams[i].lj.c6;
					complex.lj_params.push_back(lj_params);
					break;
				}//end switch
			}//end loading interaction data.

			//load molecule data
			for(size_t mol_block = 0;mol_block <mtop.nmolblock;mol_block++)
			{
				if(mtop.moltype[mtop.molblock[mol_block].type].name != 0 && (strcmp(*mtop.moltype[mtop.molblock[mol_block].type].name,"SOL") == 0/* || strcmp(*mtop.moltype[mtop.molblock[mol_block].type].name,"CL") == 0*/))// NO SOLVENTS!!!
					continue;
				for(size_t ith_mol_block = 0;ith_mol_block<mtop.molblock[mol_block].nmol;ith_mol_block++)
				{
					gmx_moltype_t& mol = mtop.moltype[mtop.molblock[mol_block].type];
					forcefield_t * curr_field;
					std::vector<mmpbsa::atom_t>* curr_atom_list;
					if(mol_block == 0)
					{
						curr_field = &split_ff[0][MMPBSAState::RECEPTOR];
						curr_atom_list = &atom_lists[0][MMPBSAState::RECEPTOR];
					}
					else
					{
						curr_field = &split_ff[0][MMPBSAState::LIGAND];
						curr_atom_list = &atom_lists[0][MMPBSAState::LIGAND];
					}

					//load atom data
					for(size_t atom_idx = 0;atom_idx<mol.atoms.nr;atom_idx++)
					{
						t_atom& atom = mol.atoms.atom[atom_idx];
						mmpbsa::atom_t new_atom;
						new_atom.atom_type = atom.typeB;
						new_atom.atomic_number = atom.atomnumber;
						new_atom.charge = charge_units*atom.qB;
						new_atom.name = (*mol.atoms.atomname[atom_idx]);
						for(size_t ex_idx = mol.excls.index[atom_idx];ex_idx < mol.excls.index[atom_idx+1];ex_idx++)
							new_atom.exclusion_list.insert(mol.excls.a[ex_idx]);
						curr_atom_list->push_back(new_atom);

						//exclusion list may need to be revised if there is an atom index offset for the complex
						if(atom_offset != 0)
						{
							new_atom.exclusion_list.clear();
							for(size_t ex_idx = mol.excls.index[atom_idx];ex_idx < mol.excls.index[atom_idx+1];ex_idx++)
								new_atom.exclusion_list.insert(mol.excls.a[ex_idx] + atom_offset);
						}
						atom_lists[0][MMPBSAState::COMPLEX].push_back(new_atom);
					}
					//load interaction indices
					for(size_t erg_type = 0;erg_type < F_NRE;erg_type++)
					{
						for(size_t erg_idx = 0;erg_idx < mol.ilist[erg_type].nr + mol.ilist[erg_type].nr_nonperturbed;erg_idx += interaction_function[erg_type].nratoms + 1)
						{
							switch(erg_type)
							{
							case F_BONDS:
							{
								bond_t new_bond;
								new_bond.atom_i = mol.ilist[erg_type].iatoms[erg_idx + 1];
								new_bond.atom_j = mol.ilist[erg_type].iatoms[erg_idx + 2];
								new_bond.bond_energy = &complex.bond_energy_data[mol.ilist[erg_type].iatoms[erg_idx] - function_start_index[erg_type]];
								curr_field->bonds_with_H.push_back(new_bond);

								//add index offset if needs for complex
								new_bond.atom_i += atom_offset;
								new_bond.atom_j += atom_offset;
								complex.bonds_with_H.push_back(new_bond);
								break;
							}
							case F_ANGLES:
							{
								angle_t new_bond;
								new_bond.atom_i = mol.ilist[erg_type].iatoms[erg_idx + 1];
								new_bond.atom_j = mol.ilist[erg_type].iatoms[erg_idx + 2];
								new_bond.atom_k = mol.ilist[erg_type].iatoms[erg_idx + 3];
								new_bond.angle_energy = &complex.angle_energy_data[mol.ilist[erg_type].iatoms[erg_idx] - function_start_index[erg_type]];
								curr_field->angles_with_H.push_back(new_bond);

								//add index offset if needs for complex
								new_bond.atom_i += atom_offset;
								new_bond.atom_j += atom_offset;
								new_bond.atom_k += atom_offset;
								complex.angles_with_H.push_back(new_bond);
								break;
							}
							case F_PDIHS: case F_PIDIHS:
							{
								dihedral_t new_dihedral;
								new_dihedral.atom_i = mol.ilist[erg_type].iatoms[erg_idx + 1];
								new_dihedral.atom_j = mol.ilist[erg_type].iatoms[erg_idx + 2];
								new_dihedral.atom_k = mol.ilist[erg_type].iatoms[erg_idx + 3];
								new_dihedral.atom_l = mol.ilist[erg_type].iatoms[erg_idx + 4];
								new_dihedral.lj.c12 = new_dihedral.lj.c6 = 0;
								new_dihedral.dihedral_energy = &complex.dihedral_energy_data[mol.ilist[erg_type].iatoms[erg_idx] - function_start_index[erg_type]];
								curr_field->dihedrals_with_H.push_back(new_dihedral);

								new_dihedral.atom_i += atom_offset;
								new_dihedral.atom_j += atom_offset;
								new_dihedral.atom_k += atom_offset;
								new_dihedral.atom_l += atom_offset;
								complex.dihedrals_with_H.push_back(new_dihedral);
								break;
							}
							case F_LJ14:
								dihedral_t new_dihedral;
								new_dihedral.atom_i = mol.ilist[erg_type].iatoms[erg_idx + 1];
								new_dihedral.atom_l = mol.ilist[erg_type].iatoms[erg_idx + 2];
								new_dihedral.atom_j = new_dihedral.atom_k = -1;
								new_dihedral.lj.c12 = nm2angst_12*joules2cal*mtop.ffparams.iparams[mol.ilist[erg_type].iatoms[erg_idx]].lj14.c12A;
								new_dihedral.lj.c6 = nm2angst_6*joules2cal*mtop.ffparams.iparams[mol.ilist[erg_type].iatoms[erg_idx]].lj14.c6A;
								new_dihedral.dihedral_energy = 0;
								curr_field->dihedrals_with_H.push_back(new_dihedral);

								new_dihedral.atom_i += atom_offset;
								new_dihedral.atom_l += atom_offset;
								complex.dihedrals_with_H.push_back(new_dihedral);
								break;
							case F_LJ:
								break;
							default:
								std::cerr << "get_gromacs_forcefield: Warning: Non supported gromacs interaction type: " << erg_type << std::endl;
								break;
							}//end switch
						}//end energy loop
					}//end outer energy loop
					atom_offset += mol.atoms.nr;
				}
			}
		}
	}

	//done_state(&state);
	sfree(f);

	//create mol_list
	size_t complex_size = atom_lists[0][MMPBSAState::COMPLEX].size();
	mol_list.resize(complex_size,MMPBSAState::LIGAND);
	mol_list[std::slice(0,atom_lists[0][MMPBSAState::RECEPTOR].size(),1)] = MMPBSAState::RECEPTOR;

	//All molecules must share the LJ paramters
	receptor.lj_params = ligand.lj_params = complex.lj_params;
}

#if 0 // crap tpr stuff


mmpbsa::SanderParm* mmpbsa_io::gmxtpr2parmtop(const std::string& filename)
{
	return gmxtpr2parmtop(filename.c_str());
}

mmpbsa::SanderParm* mmpbsa_io::gmxtpr2parmtop(std::iostream& gmxtop)
{
	throw mmpbsa::MMPBSAException("mmpbsa_io::gmxtop2parmtop: Not Implemented.");
}

std::valarray<size_t>& array_gmx2sander(mmpbsa::SanderParm* sp,const size_t& gmx_bond_type)
{
	if(sp == 0)
		throw mmpbsa::MMPBSAException("array_gmx2sander: Null pointer given for SanderParm.",mmpbsa::NULL_POINTER);

	switch(gmx_bond_type)
	{
	case F_BONDS: case F_G96BONDS:
		return sp->bonds_inc_hydrogen;
	case F_ANGLES: case F_G96ANGLES:
		return sp->angles_inc_hydrogen;
	case F_PDIHS: /*case F_IDIHS:*/
		return sp->dihedrals_inc_hydrogen;
	case F_LJ14:
		throw mmpbsa::MMPBSAException("array_gmx2sander: LJ1-4 interactions use multiple arrays.");
	default:
	{
		std::ostringstream error;
		error << "array_gmx2sander: Unsupported gromacs interaction type: " << gmx_bond_type
				<< std::endl << "See include/gromacs/types/idef.h for list of interaction types.";
		throw mmpbsa::MMPBSAException(error);
	}
	}

}

void insert_gromacs_LJ14(mmpbsa::SanderParm* sp,const t_ilist* ilist,const size_t& i, const mmpbsa_io::gromacs_idx_offsets& offsets)
{
	size_t type,a_i,a_j;

	if(ilist == 0)
		throw mmpbsa::MMPBSAException("insert_gromacs_LJ14: List of interaction parameters is needed, but a null pointer was provided.",mmpbsa::NULL_POINTER);

	type = ilist->iatoms[i*3];
	a_i = ilist->iatoms[i*3+1];
	a_j = ilist->iatoms[i*3+2];

	sp->nonbonded_parm_indices[a_i*sp->natom + a_j] = type - offsets.f_lj14 + 1/* amber adds one */;
	sp->nonbonded_parm_indices[a_j*sp->natom + a_i] = type - offsets.f_lj14 + 1/* amber adds one */;
	sp->nonbonded_parm_mask[a_j*sp->natom + a_i] = false;
	sp->nonbonded_parm_mask[a_i*sp->natom + a_j] = false;
}

void insert_molecule(mmpbsa::SanderParm* sp,const gmx_molblock_t* molblock, const gmx_moltype_t* moltype, const gmx_ffparams_t& ffparams,
		mmpbsa_io::gromacs_idx_offsets& offsets)
{
	const t_atoms* atom_list = &moltype->atoms;
	char*** atomnames = atom_list->atomname;
	size_t atom_0 = offsets.atom,funct_natoms,arr_beginning,exclusion_counter = 0;
	std::vector<size_t> bond_types = mmpbsa_io::allowed_gmx_energies();
	std::vector<size_t>::const_iterator bond_type;
	const t_ilist* ilist;
	int last_residue = -1;

	//atom data
	for(size_t i = 0;i<atom_list->nr;i++)
	{
		const t_atom& atom = atom_list->atom[i];
		sp->atom_names[atom_0 + i] = *(atomnames[i]);
		sp->charges[atom_0 + i] = atom.qB;
		sp->masses[atom_0 + i] = atom.mB;
		//add lennard jones coefficient indices later //sp->atom_type_indices[atom_0 + i] = atom.typeB;
		if(last_residue != atom.resind)
		{
			sp->residue_pointers[offsets.residue + atom.resind] = atom_0+ i + 1;//Sander uses a 1-indexed array (will be 0-indexed later by EmpEnerFun)
			sp->residue_labels[atom.resind+offsets.residue] = (atom_list->resinfo[atom.resind].name != 0) ? *atom_list->resinfo[atom.resind].name : "DRC";
			last_residue = atom.resind;
		}
		sp->iptres = offsets.residue + atom.resind;//this will guarantee that iptres is the last solute residue, since solvents will not be loaded.

		//update exclusion list
		for(size_t j = moltype->excls.index[i],k=0; j<moltype->excls.index[i+1]-1;j++)
		{
			sp->excluded_atoms_list[exclusion_counter + k++] = moltype->excls.a[j];
			sp->number_excluded_atoms[atom_0 + i] = k;
		}
		exclusion_counter += sp->number_excluded_atoms[atom_0 + i];

	}

	//interaction data
	//SanderParm data arrays are of the form
	//array = {i-th atom, j-th atom,...,force index}
	//whereas gromacs data arrays are of the form
	//array = {force index, i-th atom, j-th atom,...}
	//This loop will find the sanderparm array corresponding to each gromacs bond type
	//and insert values in the correct order.
	for(bond_type = bond_types.begin();bond_type != bond_types.end();bond_type++)
	{
		ilist = &moltype->ilist[*bond_type];
		if(ilist->iatoms == 0)
			continue;
		funct_natoms = interaction_function[*bond_type].nratoms;
		std::valarray<size_t>* interaction_array;
		size_t& funct_offset = mmpbsa_io::get_gmxfunct_offset(offsets,*bond_type);
		size_t& array_offset = mmpbsa_io::get_gmxarray_offset(offsets,*bond_type);

		try{
			interaction_array = &(array_gmx2sander(sp,*bond_type));
		}
		catch(mmpbsa::MMPBSAException mmpbsae)
		{
			if(*bond_type != F_LJ14)
				throw mmpbsae;
			interaction_array = 0;
		}

		//Iterate through interaction lists and update SanderParm array.
		for(size_t i = 0,i_0;i<ilist->nr + ilist->nr_nonperturbed;i += funct_natoms + 1)
		{
			i_0 = i/(funct_natoms+1);//i_0 represents the first entry in the current set if ilist values.
			arr_beginning = (funct_natoms+1)*(array_offset + i_0);//beginning of the i-th atom's data in the sanderparm array
			if(interaction_array != 0 && arr_beginning + funct_natoms >= interaction_array->size())
			{
				std::ostringstream error;
				error << "insert_molecule: interaction loop exceeded array size."
						<< std::endl << "Size: " << interaction_array->size()
						<< " Requested value: " << arr_beginning + funct_natoms;
				throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);
			}

			//Insert interaction data into correct SanderParm array.
			if(*bond_type == F_LJ14)
				insert_gromacs_LJ14(sp,ilist,i_0,offsets);
			else
			{
				for(size_t j = 0;j<funct_natoms;j++)
					(*interaction_array)[arr_beginning + j] = ilist->iatoms[i_0*(funct_natoms+1) + (j+1 % (funct_natoms+1))];
				(*interaction_array)[arr_beginning + funct_natoms] = ilist->iatoms[i_0*(funct_natoms+1)]-funct_offset;//shift type value to value of index in sanderparm array.
			}
		}
		array_offset += ilist->nr + ilist->nr_nonperturbed;
	}
	offsets.atom += molblock->natoms_mol;
	offsets.residue += atom_list->nres;
}

void insert_energy_check_size(const std::valarray<mmpbsa_t>& arr,const size_t& position,const std::string& function_name = "insert_energy_check_size");
void insert_energy_check_size(const std::valarray<mmpbsa_t>& arr,const size_t& position,const std::string& function_name)
{
	if(position >= arr.size())
	{
		std::ostringstream error;
		error << function_name << ": interaction loop exceeded array size."
				<< std::endl << "Size: " << arr.size()
				<< " Requested value: " << position;
		throw mmpbsa::MMPBSAException(error,mmpbsa::DATA_FORMAT_ERROR);

	}
}

bool insert_energy(mmpbsa::forcefield_t* ff,const t_functype& functype,const t_iparams params, mmpbsa_io::gromacs_idx_offsets& offsets,const std::vector<size_t>& allowed_energies)
{
	//Is this function type allowed? If not, return false. This way, if we don't care, we can ignore the fact.
	bool is_allowed = false;
	size_t& array_pos = mmpbsa_io::get_gmxarray_offset(offsets,functype);

	switch(functype)
	{
	case F_BONDS: case F_G96BONDS:
		sp->bond_equil_values[array_pos] = params.harmonic.rA;
		insert_energy_check_size(sp->bond_force_constants,array_pos,"insert_energy");
		sp->bond_force_constants[array_pos++] = params.harmonic.krA/2;//gromacs does not put 1/2 into constant. Sander Does. (for energy equation 1/2kx^2)
		break;
	case F_ANGLES: case F_G96ANGLES:
		insert_energy_check_size(sp->angle_equil_values,array_pos,"insert_energy");
		sp->angle_equil_values[array_pos] = params.harmonic.rA*MMPBSA_DEG_TO_RAD;//gromacs angles are in degrees :-/
		insert_energy_check_size(sp->angle_force_constants,array_pos,"insert_energy");
		sp->angle_force_constants[array_pos++] = params.harmonic.krA/2;//gromacs does not put 1/2 into constant. Sander Does. (for energy equation 1/2kx^2)
		break;
	case F_PDIHS:
		insert_energy_check_size(sp->dihedral_force_constants,array_pos,"insert_energy");
		sp->dihedral_force_constants[array_pos] = params.pdihs.cpA;
		insert_energy_check_size(sp->dihedral_phases,array_pos,"insert_energy");
		sp->dihedral_phases[array_pos] = params.pdihs.phiA;
		insert_energy_check_size(sp->dihedral_periodicities,array_pos,"insert_energy");
		sp->dihedral_periodicities[array_pos++] = params.pdihs.mult;
		break;
	case F_LJ14:
		insert_energy_check_size(sp->lennard_jones_acoefs,array_pos,"insert_energy");
		sp->lennard_jones_acoefs[array_pos] = params.lj14.c6A;
		insert_energy_check_size(sp->lennard_jones_bcoefs,array_pos,"insert_energy");
		sp->lennard_jones_bcoefs[array_pos++] = params.lj14.c12A;
		break;
	default:
		return false;
	}

	return true;
}

mmpbsa::SanderParm* mmpbsa_io::gmxtpr2parmtop(const char* fn)
{
	if(true)
		throw "gmxtpr2parmtop is super depricated";

	mmpbsa::SanderParm* returnMe;
	mmpbsa_io::gromacs_idx_offsets offsets;
	size_t molecule_counter;
	FILE *gp;
	  int         fp,indent,i,j,**gcount,atot;
	  t_state     state;
	  rvec        *f=NULL;
	  t_inputrec  ir;
	  t_tpxheader tpx;
	  t_topology top;
	  gmx_mtop_t  mtop;
	  gmx_groups_t *groups;
	  gmx_ffparams_t* force_field;
	  const char* mdpfn = 0;
	  bool bSysTop,bShowNumbers;
	  bSysTop = false;
	  bShowNumbers = true;

	  if(fn == 0 || strlen(fn) == 0)
		  throw mmpbsa::MMPBSAException("mmpbsa_io::gmxtpr2parmtop: Filename is required.",mmpbsa::DATA_FORMAT_ERROR);

	  returnMe = new mmpbsa::SanderParm;
	  printf("Make sure to check DIELECTRIC CONSTANTS AND GROMACS UNITS!!!!");
	  const std::vector<size_t>& interaction_types = mmpbsa_io::allowed_gmx_energies();
	  std::vector<size_t>::const_iterator interaction_type;

//#define DO_GMXDUMP 1

#if DO_GMXDUMP//old
	  read_tpxheader(fn,&tpx,TRUE,NULL,NULL);
	  read_tpx_state(fn,
			 tpx.bIr  ? &ir : NULL,
			 &state,tpx.bF ? f : NULL,
			 tpx.bTop ? &mtop: NULL);
#endif
#ifndef DO_GMXDUMP
	  read_tpx(fn,NULL,NULL,&atot,0,0,0,&mtop);
#endif
#if 0
	  read_tpxheader(fn,&tpx,TRUE,NULL,NULL);
	  read_tpx_state(fn,NULL,NULL,NULL,&mtop);
#endif


#if 0//don't need
	  if (mdpfn && tpx.bIr) {
	    gp = gmx_fio_fopen(mdpfn,"w");
	    pr_inputrec(gp,0,NULL,&(ir),TRUE);
	    gmx_fio_fclose(gp);
	  }
#endif

	  if (!mdpfn) {
		  //top = gmx_mtop_t_to_t_topology(&mtop);
	    if (available(stdout,&tpx,0,fn)) {
#ifdef DO_GMXDUMP
	      indent=0;
	      indent=pr_title(stdout,indent,fn);
	      pr_inputrec(stdout,0,"inputrec",tpx.bIr ? &(ir) : NULL,FALSE);

	      indent = 0;
	      pr_header(stdout,indent,"header",&(tpx));
#endif
	      //setup constants
	      init(offsets);
	      returnMe->titles = (mtop.name != 0 && *mtop.name != 0) ? *mtop.name : "";
	      returnMe->natom = mtop.natoms;
	      returnMe->ntypes = returnMe->natom;//to go to amber lj interaction matrix, type will equal atom id.
	      returnMe->nres = 0;
	      for(size_t mb = 0;mb<mtop.nmolblock;mb++)
	      {
	    	  if(strcmp(*mtop.moltype[mtop.molblock[mb].type].name,"SOL") == 0)//remove/ignore solvent
	    		  returnMe->natom -= mtop.molblock[mb].nmol*mtop.molblock[mb].natoms_mol;
	    	  else
	    		  returnMe->nspm +=  mtop.molblock[mb].nmol;
	      }
	      for(size_t mt = 0;mt<mtop.nmoltype;mt++)
	      {
	    	  if(mtop.moltype[mt].name != 0 && strcmp(*mtop.moltype[mt].name,"SOL") == 0)//ignore solvent.
	    		  continue;

	    	  if(mtop.moltype[mt].excls.nra != 0)
	    		  returnMe->nnb += mtop.moltype[mt].excls.nra;
	    	  returnMe->nres += mtop.moltype[mt].atoms.nres;
	    	  for(interaction_type = interaction_types.begin();interaction_type != interaction_types.end();interaction_type++)
	    	  {
	    		  mmpbsa_io::get_gmxfunct_offset(offsets,*interaction_type) = -1;//max value of size_t since it is unsigned.
	    		  switch(*interaction_type)
	    		  	{
	    		  	case F_BONDS: case F_G96BONDS:
	    		  		returnMe->nbonh += mtop.moltype[mt].ilist[*interaction_type].nr/3;
	    		  		break;
	    		  	case F_ANGLES: case F_G96ANGLES:
	    		  		returnMe->ntheth += mtop.moltype[mt].ilist[*interaction_type].nr/4;
	    		  		break;
	    		  	case F_PDIHS:/* case F_IDIHS:*/
	    		  		returnMe->nphih += mtop.moltype[mt].ilist[*interaction_type].nr/5;
	    		  		break;
	    		  	case F_LJ14:
	    		  		break;
	    		  	default:
	    		  	{
	    		  		std::ostringstream error;
	    		  		error << "mmpbsa_io::gmxtpr2parmtop: Unsupported gromacs interaction type: " << *interaction_type
	    		  				<< std::endl << "See include/gromacs/types/idef.h for list of interaction types.";
	    		  		throw mmpbsa::MMPBSAException(error);
	    		  	}
	    		  	}
	    	  }

	      }
	      for(size_t i = 0;i<mtop.ffparams.ntypes;i++)
	      {
	    	  switch(mtop.ffparams.functype[i])
	    	  {
	    	  case F_LJ14/* use this to set ntypes??? */:
	    		  break;
	    	  case F_BONDS:case F_G96BONDS:
	    		  returnMe->numbnd++;
	    		  break;
	    	  case F_ANGLES: case F_G96ANGLES:
	    		  returnMe->numang++;
	    		  break;
	    	  case F_PDIHS:/* case F_IDIHS: */
	    		  returnMe->nptra++;
	    		  break;
	    	  default:
	    		  continue;
	    	  }
	    	  if(i < mmpbsa_io::get_gmxfunct_offset(offsets,mtop.ffparams.functype[i]))//update offsets
	    		  mmpbsa_io::get_gmxfunct_offset(offsets,mtop.ffparams.functype[i]) = i;
	      }

	      returnMe->initialize_arrays();//initializes arrays based on size data in sanderparm.

	      force_field = &mtop.ffparams;//when loading force constants divide all by 2 except for proper dihedral
	      molecule_counter = 0;

	      //setup LJ atom type indices. Gromacs does not use the same types for LJ 1-4 calculations and sander.
	      //So here the "atom type" is simply the atom's index. This will make the array's bigger, but at least it works.
	      //Will use more memory, but not more time.
	      for(size_t i = 0;i<returnMe->natom;i++)
	    	  returnMe->nonbonded_parm_indices[i] = i+1;

	      //Insert molecules into SanderParm
	      for(size_t mb = 0;mb<mtop.nmolblock;mb++)
	    	  for(size_t j = 0;j<mtop.molblock[mb].nmol;j++)//normally will only run one j-loop. Solvents are excluded.
	    	  {
	    		  if(strcmp(*mtop.moltype[mtop.molblock[mb].type].name,"SOL") == 0)
	    			  continue;
	    		  insert_molecule(returnMe,&mtop.molblock[mb],&mtop.moltype[mb],*force_field, offsets);
	    		  returnMe->atoms_per_molecule[molecule_counter++] = mtop.molblock[mb].natoms_mol;
	    	  }

	      //Insert Force Field data into SanderParm
	      init(offsets);//offsets will keep track of position in arrays
	      for(size_t ff = 0;ff<mtop.ffparams.ntypes;ff++)
	      {
	    	  insert_energy(returnMe,mtop.ffparams.functype[ff],mtop.ffparams.iparams[ff],offsets, interaction_types);
	      }

#if DO_GMXDUMP //don't need

	      if (!bSysTop)
	    	  pr_mtop(stdout,indent,"topology",&(mtop),bShowNumbers);
	      else
	    	  pr_top(stdout,indent,"topology",&(top),bShowNumbers);
	      pr_rvecs(stdout,indent,"box",tpx.bBox ? state.box : NULL,DIM);
	      pr_rvecs(stdout,indent,"box_rel",tpx.bBox ? state.box_rel : NULL,DIM);
	      pr_rvecs(stdout,indent,"boxv",tpx.bBox ? state.boxv : NULL,DIM);
	      pr_rvecs(stdout,indent,"pres_prev",tpx.bBox ? state.pres_prev : NULL,DIM);
	      pr_rvecs(stdout,indent,"svir_prev",tpx.bBox ? state.svir_prev : NULL,DIM);
	      pr_rvecs(stdout,indent,"fvir_prev",tpx.bBox ? state.fvir_prev : NULL,DIM);
	      /* leave nosehoover_xi in for now to match the tpr version */
	      pr_doubles(stdout,indent,"nosehoover_xi",state.nosehoover_xi,state.ngtc);
	      /*pr_doubles(stdout,indent,"nosehoover_vxi",state.nosehoover_vxi,state.ngtc);*/
	      /*pr_doubles(stdout,indent,"therm_integral",state.therm_integral,state.ngtc);*/
	      pr_rvecs(stdout,indent,"x",tpx.bX ? state.x : NULL,state.natoms);
	      pr_rvecs(stdout,indent,"v",tpx.bV ? state.v : NULL,state.natoms);
	      if (state,tpx.bF) {
		pr_rvecs(stdout,indent,"f",f,state.natoms);
	      }
#endif
	    }

#if DO_GMXDUMP//don't need
	    groups = &mtop.groups;

	    snew(gcount,egcNR);
	    for(i=0; (i<egcNR); i++)
	      snew(gcount[i],groups->grps[i].nr);

	    for(i=0; (i<mtop.natoms); i++) {
	      for(j=0; (j<egcNR); j++)
		gcount[j][ggrpnr(groups,j,i)]++;
	    }
	    printf("Group statistics\n");
	    for(i=0; (i<egcNR); i++) {
	      atot=0;
	      printf("%-12s: ",gtypes[i]);
	      for(j=0; (j<groups->grps[i].nr); j++) {
		printf("  %5d",gcount[i][j]);
		atot+=gcount[i][j];
	      }
	      printf("  (total %d atoms)\n",atot);
	      sfree(gcount[i]);
	    }
	    sfree(gcount);
#endif

	  }
	  done_state(&state);
	  sfree(f);
	  return returnMe;
}
#endif//depricated


