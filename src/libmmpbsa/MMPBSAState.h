/**
 * @class mmpbsa::MMPBSAState
 * @brief State structure for Molecular Dynamics and MMPBSA
 *
 * Saves the state of the MD or MMPBSA calculation, which might
 * be needed to restart calculations.
 */
#ifndef MMPBSASTATE_H
#define MMPBSASTATE_H

#include <map>
#include <string>
#include <vector>

#include "libmmpbsa/SanderInterface.h"
#include "libmmpbsa/MeadInterface.h"

#define SANDER_MDOUT_TYPE "mdout"
#define SANDER_PRMTOP_TYPE "prmtop"
#define MMPBSA_OUT_TYPE "mmpbsa_out"
#define RADII_TYPE "radii"
#define SANDER_INPCRD_TYPE "inpcrd"
#define CHECKPOINT_FILE_TYPE "checkpoint"

//These tags are used to indicate the current state of the mdmmpbsa program
#define MMPBSASTATE_TAG "mmpbsa_state"
#define MDSTATE_TAG "moldyn_state"

namespace mmpbsa{

class MMPBSAState
{
public:
    std::vector<size_t> receptorStartPos;///<Starting positions for receptors. End positions are deduced from the parmtop file.
    std::vector<size_t> ligandStartPos;///<Starting positions for ligands. End positions are deduced from the parmtop file.
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
    enum SimProcess{MMPBSA,SANDER} currentProcess;
    int placeInQueue;//zero ordered place in queue.
    float weight;//factor by which to multiple time in queue(i.e. weight = 1 means process will take equal time as other processes).

    bool trustPrmtop;///<Flag to indicate if the sanity check of the SanderParm object should be ignored. This is not suggested, but if the sanity check fails and one *does* believe it should work, this is provided as a work around, for whatever reason might arise.
    bool keep_traj_in_mem;///<Flag to indicate whether or not the trajectory stream should stay in memory. Default: false

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

bool has_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state);
const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state) throw (mmpbsa::MMPBSAException);


#endif//MMPBSASTATE_H

