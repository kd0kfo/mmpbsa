#include "GromacsReader.h"

void mmpbsa_io::load_gmx_trr(const std::string& filename,std::valarray<mmpbsa_t>& crds,size_t frame_number)
{
  t_fileio    *fpread/* ,*fpwrite */;
  size_t         nframe/*,indent*/;
  rvec        *x,*v,*f;
  matrix      box;
  t_trnheader trn;
  gmx_bool        bOK;
  size_t natoms;

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
