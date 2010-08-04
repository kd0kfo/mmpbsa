#ifndef MMPBSA_H
#define	MMPBSA_H

#include <cstdlib>
#include <iostream>
#include <valarray>
#include <fstream>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libmmpbsa/mmpbsa_exceptions.h"
#include "libmmpbsa/mmpbsa_io.h"
#include "libmmpbsa/EnergyInfo.h"
#include "libmmpbsa/EmpEnerFun.h"
#include "libmmpbsa/MeadInterface.h"
#include "libmmpbsa/XMLParser.h"
#include "libxml/tree.h"//include libxml2 in Includes path.
#include "libmmpbsa/SanderInterface.h"

#include "MEAD/FinDiffMethod.h"

#ifdef __USE_BOINC__
#include "boinc_api.h"
#include "str_util.h"
#else
#define EXIT_CHILD_FAILED 195
#endif

double timeAtPreviousCheckpoint;

class MMPBSAState
{
public:
    std::vector<size_t> receptorStartPos;///<Starting positions for receptors. End positions are deduced from the parmtop file.
    std::vector<size_t> ligandStartPos;///<Starting positions for ligands. End positions are deduced from the parmtop file.
    std::vector<size_t> snapList;///<List of snapshots to be used in the calculations. If the list is empty, all snapshots are used.

    std::string outputFilename;
    std::string prmtopFilename;
    std::string trajFilename;
    std::string radiiFilename;
    std::string checkpointFilename;

    size_t checkpointCounter;
    size_t currentSnap;

    double fractionDone;

    bool MDOnly;

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
    int placeInQueue;//zero orderd place in queue.
    float weight;//factor by which to multiple time in queue(i.e. weight = 1 means process will take equal time as other processes).

    bool trustPrmtop;///<Flag to indicate if the sanity check of the SanderParm object should be ignored. This is not suggested, but if the sanity check fails and one *does* believe it should work, this is provided as a work around, for whatever reason might arise.
    /**
     * Stores variables needed to restart the program. This is needed for running
     * on systems where the program may be closed and restarted again, i.e. BOINC Grid
     */
    MMPBSAState();
    MMPBSAState(const MMPBSAState& orig);
    virtual ~MMPBSAState(){}

    MMPBSAState& operator=(const MMPBSAState& orig);
};

std::vector<MMPBSAState> processQueue;

/**
 * Pulls everything together to perform the MMPBSA calculations.
 * Calls the parseArgs to load information from the command line.
 * 
 * @param argc
 * @param argv
 * @return
 */
int mmpbsa_run(MMPBSAState& currState, mmpbsa::MeadInterface& mi);

int sander_run(MMPBSAState& currState,mmpbsa::SanderInterface& si);

void updateMMPBSAProgress(MMPBSAState& currState,const double& increment);

void printSnapshot(const mmpbsa::EMap& complexEMap, const mmpbsa::EMap& receptorEMap,
        const mmpbsa::EMap& ligandEMap, std::fstream& outFile);

/**
 * Takes the command line arguments and performs the necessary functions.
 * 
 * @param argc
 * @param argv
 * @param mi
 */
std::map<std::string,std::string> parseArgs(int argc, char** argv);

/**
 * When a command line argument provides data or information(to the right of an
 * equal sign), this method decides what to do with it.
 *
 * @param arg
 * @param mi
 */
int parseParameter(std::map<std::string,std::string> args, MMPBSAState& currState, mmpbsa::MeadInterface& mi);
int parseParameter(std::map<std::string,std::string> args, MMPBSAState& currState, mmpbsa::SanderInterface& si);

/**
 * When a command line argument toggles a program flag, this method makes the
 * proper changes to the program.
 * 
 * @param flag
 * @param mi
 */
int parseFlag(std::string flag, MMPBSAState& currState, mmpbsa::MeadInterface& mi);
int parseFlag(std::string flag, MMPBSAState& currState, mmpbsa::SanderInterface& si);

/**
 * When a parameter has a list of variables (separated by commas), this method
 * tokenizes that list and loades it into the array, in left to right order.
 * 
 * @param values
 * @param array
 */
int loadListArg(const std::string& values,std::vector<size_t>& array, const size_t& offset = 0);

std::vector<MMPBSAState> getQueueFile(int argc,char** argv);

/**
 * Creates a string about the usage of the program, listing the parameters and
 * flags. This string should be sent to STDOUT, if "--help" is provided as a flag.
 * 
 * @return
 */
std::string helpString();

/**
 * Runs the boinc initialization process, if compiled with BOINC. Otherwise,
 * nothing is done.
 * 
 * @return
 */
int mmpbsa_boinc_init();

/**
 * Check with the BOINC client to see if anything should be done, if compiled
 * with BOINC. Otherwise, nothing is done.
 * 
 * @param si
 */
void poll_boinc_messages(mmpbsa::SanderInterface& si);



/**
 * Saves the state of the MMPBSA Program. If BOINC is being used, BOINC checkpointing
 * is performed.
 */
void checkpoint_out(MMPBSAState& saveState, mmpbsa_utils::XMLParser& xmlDoc);
void checkpoint_mmpbsa(MMPBSAState& saveState);
void checkpoint_sander(MMPBSAState& saveState, mmpbsa::SanderInterface& si);
\
/**
 * Load the last MMPBSA state. Returns true is all of the parameters in the
 * checkpoint file were used.
 *
 */
bool restart_mmpbsa(MMPBSAState& restartState);
bool restart_sander(MMPBSAState& restartState, mmpbsa::SanderInterface& si);

/**
 * Generates a sample queue file and writes it to the provided file. The sample
 * file is an XML document similar to that which the program takes and includes
 * parameters used by both molecular dynamics and MMPBSA.
 * 
 * @param filename
 */
void sampleQueue(const std::string& filename);

/**
 * IF BOINC is being used, the progress is reported to the client based on CPU
 * and precentages stored in the MMPBSAState's
 */
void report_boinc_progress();


#endif	/* MMPBSA_H */

