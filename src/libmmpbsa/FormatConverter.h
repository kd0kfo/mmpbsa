#ifndef FORMATCONVERTER_H
#define FORMATCONVERTER_H

#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"
#include "SanderParm.h"

namespace mmpbsa_io{

typedef struct {
	size_t curr_snap;
	//for sander trajectories
	std::iostream* sander_crd_stream;
	mmpbsa::SanderParm* sander_parm;

	//for gromacs trajectories
	std::string* gromacs_filename;
}trajectory_t;

bool get_next_snap(mmpbsa_io::trajectory_t& traj, std::valarray<mmpbsa_t>& snapshot);

void seek(mmpbsa_io::trajectory_t& traj,const size_t& snap_pos);

void default_trajectory(mmpbsa_io::trajectory_t& traj);
mmpbsa_io::trajectory_t open_trajectory(const std::string& filename);
void destroy_trajectory(mmpbsa_io::trajectory_t& traj);

bool eof(trajectory_t& traj);

std::string get_traj_title(mmpbsa_io::trajectory_t& traj);

mmpbsa::SanderParm* gmxtop2parmtop(const std::string& filename);
mmpbsa::SanderParm* gmxtop2parmtop(std::iostream& gmxtop);

}//end namespace mmpbsa_io


#endif // TRAJECTORYREADER_H

