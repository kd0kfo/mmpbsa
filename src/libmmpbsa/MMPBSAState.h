/**
 * @class mmpbsa::MMPBSAState
 * @brief State structure for Molecular Dynamics and MMPBSA
 *
 * Saves the state of the MD or MMPBSA calculation, which might
 * be needed to restart calculations.
 *
 * MMPBSAState stores all of the filenames for the program.
 * This is done using a map of file types (c.f. globals.h)
 * and file names.
 */
#ifndef MMPBSASTATE_H
#define MMPBSASTATE_H

#include <map>
#include <string>
#include <vector>

#include "globals.h"
#include "SanderInterface.h"
#include "MeadInterface.h"


namespace mmpbsa{
//forward declarations

class MMPBSAState
{
public:
    std::set<size_t> receptorStartPos;///<Starting positions for receptors. End positions are deduced from the parmtop file.
    std::set<size_t> ligandStartPos;///<Starting positions for ligands. End positions are deduced from the parmtop file.
    std::vector<size_t> snapList;///<List of snapshots to be used in the calculations. If the list is empty, all snapshots are used.

    std::map<std::string,std::string> filename_map;///<Maps a file type flag (cf Amber Manual, Sander Section) to the filename

    size_t checkpointCounter;///<Counter of the number of times checkpoint is done. For purely informative purposes.
    size_t currentSnap;///<Gives the current snapshot being used by mmpbsa_run

    double fractionDone;///<Progress of the program as a fraction of one.

    bool MDOnly;///<Indicator that only Molecular dynamics should be performed. Generally only needed when no XML queue file is provided.
    bool savePDB;///<Indicator of whether or not pdb files should be produced for each energy fun calculated.
    mmpbsa::SanderInterface currentSI;
    mmpbsa::MeadInterface currentMI;

    /**
     * List parts of the complex on which to perform MMPBSA. END_OF_MOLECUES signifies the snapshot is finished.
     */
    enum MOLECULE{COMPLEX,RECEPTOR,LIGAND,END_OF_MOLECULES} currentMolecule;

    /**
     * Stage in the program.
     */
    enum SimProcess{MMPBSA,SANDER,MOLSURF} currentProcess;
    int placeInQueue;//zero ordered place in queue.
    float weight;//factor by which to multiple time in queue(i.e. weight = 1 means process will take equal time as other processes).

    bool trustPrmtop;///<Flag to indicate if the sanity check of the SanderParm object should be ignored. This is not suggested, but if the sanity check fails and one *does* believe it should work, this is provided as a work around, for whatever reason might arise.
    bool keep_traj_in_mem;///<Flag to indicate whether or not the trajectory stream should stay in memory. Default: false
    bool surface_area_only;///<Flag to indicate that only SA should be done in PBSA

    int verbose;///<Flag to indicate whether the program needs to be verbose. Added in version 0.12.5. Not fully implemented yet

    /**
     * Stores variables needed to restart the program. This is needed for running
     * on systems where the program may be closed and restarted again, i.e. BOINC Grid
     */
    MMPBSAState();
    MMPBSAState(const MMPBSAState& orig);
    virtual ~MMPBSAState(){}

    MMPBSAState& operator=(const MMPBSAState& orig);
};

}//namespace mmpbsa

/**
 * Determines whether the specified mmpbsa::MMPBSAState object has a
 * a record of specific file type.
 *
 * @see const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state) throw (mmpbsa::MMPBSAException);
 */
bool has_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state);

/**
 * Returns a reference to the filename of the specific file type. If
 * the file type is not listed in the mmpbsa::MMPBSAState object,
 * an mmpbsa::MMPBSAException is thrown. This can be prevented
 * if the object is first checked using, has_filename
 */
const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state) throw (mmpbsa::MMPBSAException);


#endif//MMPBSASTATE_H

