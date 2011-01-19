#include "GromacsReader.h"

void mmpbsa_io::load_gmx_trr(const std::string& filename,std::valarray<mmpbsa_t>& crds,size_t frame_number,const size_t* natom_limit)
{
  t_fileio    *fpread/* ,*fpwrite */;
  size_t         nframe/*,indent*/;
  rvec        *x,*v,*f;
  matrix      box;
  t_trnheader trn;
  gmx_bool        bOK;
  size_t natoms;

  if(filename.size() == 0)
	  throw mmpbsa::MMPBSAException("mmpbsa_io::load_gmx_trr: File name required.");

  fpread  = open_trn(filename.c_str(),"r");
  //fpwrite = open_tpx(NULL,"w");
  //gmx_fio_setdebug(fpwrite,TRUE);

  nframe = 0;
  natoms = 0;
  x = 0;
  while (fread_trnheader(fpread,&trn,&bOK)) {
    snew(x,trn.natoms);
    snew(v,trn.natoms);
    snew(f,trn.natoms);
    if (fread_htrn(fpread,&trn,
		   trn.box_size ? box : NULL,
		   trn.x_size   ? x : NULL,
		   trn.v_size   ? v : NULL,
		   trn.f_size   ? f : NULL)) {
#if 0//UNNECESSARY
      sprintf(buf,"%s frame %d",fn,nframe);
      indent=0;
      indent=pr_title(stdout,indent,buf);
      pr_indent(stdout,indent);
      fprintf(stdout,"natoms=%10d  step=%10d  time=%12.7e  lambda=%10g\n",
      	      trn.natoms,trn.step,trn.t,trn.lambda);
      if (trn.box_size)
    	  pr_rvecs(stdout,indent,"box",box,DIM);
      if (trn.x_size)
    	  pr_rvecs(stdout,indent,"x",x,trn.natoms);
      if (trn.v_size)
    	  pr_rvecs(stdout,indent,"v",v,trn.natoms);
      if (trn.f_size)
    	  pr_rvecs(stdout,indent,"f",f,trn.natoms);
#endif
    }
    else
    	std::cerr << "mmpbsa_io::load_gmx_trr: WARNING: Incomplete frame: nr " << nframe << ", t=" << trn.t << std::endl;

    if(frame_number == nframe)
    {
    	natoms = trn.natoms;
    	break;
    }
    else
    {
    	sfree(x);
    	x = 0;
    }
    sfree(v);
    sfree(f);
    nframe++;
  }
  if (!bOK)
    std::cerr << "mmpbsa_io::load_gmx_trr: WARNING: Incomplete frame header: nr " << nframe << ", t=" << trn.t << std::endl;

  close_trn(fpread);

  if(natom_limit != 0)
	  natoms = *natom_limit;

  if(natoms != 0 && x != 0)
  {
	  crds.resize(natoms*3,0);
	  for(size_t i = 0;i<natoms;i++)
		  for(size_t j = 0;j<3;j++)
			  crds[i*3+j] = x[i][j];
  }
  else
	  crds.resize(0,0);
  if(x != 0)
	  sfree(x);
}

std::vector<size_t> mmpbsa_io::allowed_gmx_energies()
{
	std::vector<size_t> returnMe;
	returnMe.push_back(F_BONDS);
	returnMe.push_back(F_G96BONDS);
	returnMe.push_back(F_ANGLES);
	returnMe.push_back(F_G96ANGLES);
	returnMe.push_back(F_PDIHS);
	//returnMe.push_back(F_IDIHS);
	returnMe.push_back(F_LJ14);
	return returnMe;
}

void init(mmpbsa_io::gromacs_idx_offsets& offset)
{
	offset.atom = 0;

	offset.f_bonds = offset.f_g96bonds = offset.f_angles = offset.f_g96angles = offset.f_pdihs = /*offset.f_idihs =*/ offset.f_lj14 = 0;
	offset.residue = 0;
	offset.bndh = offset.bnda = offset.thetah = offset.theta = offset.phih = offset.phia = offset.lj14 = 0 ;
}

size_t& mmpbsa_io::get_gmxarray_offset(mmpbsa_io::gromacs_idx_offsets& offsets,const size_t& gmx_bond_type)
{
	switch(gmx_bond_type)
	{
	case F_LJ14:
		return offsets.lj14;
	case F_BONDS: case F_G96BONDS:
		return offsets.bndh;
	case F_ANGLES: case F_G96ANGLES:
		return offsets.thetah;
	case F_PDIHS: /*case F_IDIHS:*/
		return offsets.phih;
	default:
	{
		std::ostringstream error;
		error << "mmpbsa_io::get_gmxarray_offset: Unsupported gromacs interaction type: " << gmx_bond_type
				<< std::endl << "See include/gromacs/types/idef.h for list of interaction types.";
		throw mmpbsa::MMPBSAException(error);
	}
	}
}

size_t& mmpbsa_io::get_gmxfunct_offset(mmpbsa_io::gromacs_idx_offsets& offsets,const size_t& gmx_bond_type)
{
	switch(gmx_bond_type)
	{
	case F_LJ14:
		return offsets.f_lj14;
	case F_BONDS:
		return offsets.f_bonds;
	case F_G96BONDS:
		return offsets.f_g96bonds;
	case F_ANGLES:
		return offsets.f_angles;
	case F_G96ANGLES:
		return offsets.f_g96angles;
	case F_PDIHS:
		return offsets.f_pdihs;
	/*case F_IDIHS:
		return offsets.f_idihs;*/
	default:
	{
		std::ostringstream error;
		error << "mmpbsa_io::get_gmxarray_offset: Unsupported gromacs interaction type: " << gmx_bond_type
				<< std::endl << "See include/gromacs/types/idef.h for list of interaction types.";
		throw mmpbsa::MMPBSAException(error);
	}
	}
}



