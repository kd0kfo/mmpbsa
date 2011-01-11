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

mmpbsa::SanderParm* gmxtop2parmtop(const std::string& filename)
{
	std::fstream the_file(filename,std::ios::in);
	mmpbsa::SanderParm* returnMe = gmxtop2parmtop(the_file);
	the_file.close();
	return returnMe;
}

mmpbsa::SanderParm* gmxtop2parmtop(std::iostream& gmxtop)
{

}



