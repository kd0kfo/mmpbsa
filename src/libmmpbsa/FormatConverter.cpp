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

	if(traj.sander_parm == 0 || traj.sander_crd_stream == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::get_next_snap: Sander parameters and/or sander coordinate stream is missing.",mmpbsa::DATA_FORMAT_ERROR);
	return get_next_snap(*traj.sander_crd_stream,snapshot,traj.sander_parm->natom,(traj.sander_parm->ifbox > 0));
}


void mmpbsa_io::seek(mmpbsa_io::trajectory_t& traj,const size_t& snap_pos)
{
	size_t i = traj.curr_snap;
	if(traj.sander_crd_stream != 0)
	{
		if(traj.sander_parm == 0)
			throw mmpbsa::MMPBSAException("mmpbsa_io::seek: Trajectory cannot be read without paramters. However, sander paramter object is a null pointer.",mmpbsa::NULL_POINTER);
		bool isPeriodic = (traj.sander_parm->ifbox > 0);//Are periodic boundary conditions used?
		try
		{
			for(;i<snap_pos;i++)
			{
				mmpbsa_io::skip_next_snap(*traj.sander_crd_stream,traj.sander_parm->natom,isPeriodic);

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
	traj.sander_parm = 0;

	traj.gromacs_filename = 0;
	traj.curr_snap = 0;
}

void mmpbsa_io::destroy_trajectory(mmpbsa_io::trajectory_t& traj)
{
	delete traj.sander_crd_stream;
	//don't delete sander_parm!
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
		return false;
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

bool insert_energy(mmpbsa::SanderParm* sp,const t_functype& functype,const t_iparams params, mmpbsa_io::gromacs_idx_offsets& offsets,const std::vector<size_t>& allowed_energies)
{
	//Is this function type allowed? If not, return false. This way, if we don't care, we can ignore the fact.
	bool is_allowed = false;
	for(std::vector<size_t>::const_iterator energy = allowed_energies.begin();energy != allowed_energies.end();energy++)
		if(*energy == functype)
			is_allowed = true;
	if(!is_allowed)
		return false;

	size_t& array_pos = mmpbsa_io::get_gmxarray_offset(offsets,functype);

	switch(functype)
	{
	case F_BONDS: case F_G96BONDS:
		insert_energy_check_size(sp->bond_equil_values,array_pos,"insert_energy");
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
	}

	return true;
}

mmpbsa::SanderParm* mmpbsa_io::gmxtpr2parmtop(const char* fn)
{
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



