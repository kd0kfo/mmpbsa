
#include "libmmpbsa/MMPBSAState.h"

mmpbsa::MMPBSAState::MMPBSAState()
{
    trustPrmtop = false;
    currentSnap = 0;
    currentProcess = MMPBSAState::SANDER;
    currentMolecule = MMPBSAState::COMPLEX;
    checkpointCounter = 0;
    this->fractionDone = 0;
    placeInQueue = 0;
    this->weight = 1;
    MDOnly = false;
    savePDB = false;
    keep_traj_in_mem = false;
    surface_area_only = false;
    verbose = 0;
    overwrite = false;
}

mmpbsa::MMPBSAState::MMPBSAState(const mmpbsa::MMPBSAState& orig)
{
    receptorStartPos = orig.receptorStartPos;
    ligandStartPos = orig.ligandStartPos;
    snapList = orig.snapList;
    checkpointCounter = orig.checkpointCounter;
    currentSnap = orig.currentSnap;
    fractionDone = orig.fractionDone;
    MDOnly = orig.MDOnly;
    savePDB = orig.savePDB;
    keep_traj_in_mem = orig.keep_traj_in_mem;
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
    filename_map = orig.filename_map;
    surface_area_only = orig.surface_area_only;
    verbose = orig.verbose;
    overwrite = orig.overwrite;

}

mmpbsa::MMPBSAState& mmpbsa::MMPBSAState::operator=(const mmpbsa::MMPBSAState& orig)
{
    if(this == &orig)
        return *this;

    receptorStartPos = orig.receptorStartPos;
    ligandStartPos = orig.ligandStartPos;
    snapList = orig.snapList;
    checkpointCounter = orig.checkpointCounter;
    currentSnap = orig.currentSnap;
    fractionDone = orig.fractionDone;
    MDOnly = orig.MDOnly;
    savePDB = orig.savePDB;
    keep_traj_in_mem = orig.keep_traj_in_mem;
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
    filename_map = orig.filename_map;
    surface_area_only = orig.surface_area_only;
    verbose = orig.verbose;
    overwrite = orig.overwrite;


    return *this;
}


bool has_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state)
{
	return the_state.filename_map.find(filetype) != the_state.filename_map.end();
}

const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state) throw (mmpbsa::MMPBSAException)
{
	if(the_state.filename_map.find(filetype) == the_state.filename_map.end())
		throw mmpbsa::MMPBSAException("get_filename: No filename provided for the file type " + filetype,mmpbsa::COMMAND_LINE_ERROR);
	return the_state.filename_map.find(filetype)->second;
}
