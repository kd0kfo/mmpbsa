#ifndef MMPBSA_H
#define	MMPBSA_H

#include <cstdlib>
#include <iostream>
#include <valarray>
#include <fstream>
#include <sstream>
#include <time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#define PACKAGE_STRING "MD/MMPBSA"
#endif

#include "libmmpbsa/mmpbsa_exceptions.h"
#include "libmmpbsa/mmpbsa_io.h"
#include "libmmpbsa/EnergyInfo.h"
#include "libmmpbsa/EmpEnerFun.h"
#include "libmmpbsa/EMap.h"
#include "libmmpbsa/MeadInterface.h"
#include "libmmpbsa/XMLParser.h"
#include "libmmpbsa/XMLNode.h"
#include "libmmpbsa/SanderInterface.h"
#include "libmmpbsa/MMPBSAState.h"

#include "MEAD/FinDiffMethod.h"

#ifdef USE_BOINC
#include "boinc/boinc_api.h"
#include "boinc/str_util.h"

#ifdef USE_GRAPHICS
#include "mmpbsa_graphics.h"
#include "boinc/graphics2.h"
MMPBSA_SHMEM* gshmem;
#endif

#else
#define EXIT_CHILD_FAILED 195
#endif

double timeAtPreviousCheckpoint;
double netFractionDone;
double netCPUTime;
static int mmpbsa_verbosity;


std::vector<mmpbsa::MMPBSAState> processQueue;///<Array of calculations to be run by the program.

/**
 * Pulls everything together to perform the MMPBSA calculations.
 * Calls the parseArgs to load information from the command line.
 * 
 * @param argc
 * @param argv
 * @return
 */
int mmpbsa_run(mmpbsa::MMPBSAState& currState, mmpbsa::MeadInterface& mi);

/**
 * Forks and monitors a Sander Process. Parameters are stored in the
 * mmpbsa::SanderInterface object and progress is reported by the MMPBSAState
 * object.
 *
 * If BOINC is being used, sander_run listens for messages (eg Quit or Suspend)
 *  from the BOINC client and acts accordingly.
 *
 * @param currState
 * @param si
 * @return
 */
int sander_run(mmpbsa::MMPBSAState& currState,mmpbsa::SanderInterface& si);

/**
 * Updates the progress measurement (fractionDone) for the MMPBSAState object
 * associated with the running MMPBSA calculation.
 * 
 * @param currState
 * @param increment
 */
void updateMMPBSAProgress(mmpbsa::MMPBSAState& currState,const double& increment);

/**
 * Prints the energy data for a given snapshot, using the overloaded "<<" operator
 * in EMap
 * 
 * @param complexEMap
 * @param receptorEMap
 * @param ligandEMap
 * @param outFile
 */
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
int parseParameter(std::map<std::string,std::string> args, mmpbsa::MMPBSAState& currState, mmpbsa::MeadInterface& mi);
int parseParameter(std::map<std::string,std::string> args, mmpbsa::MMPBSAState& currState, mmpbsa::SanderInterface& si);

/**
 * When a command line argument toggles a program flag, this method makes the
 * proper changes to the program.
 * 
 * @param flag
 * @param mi
 */
int parseFlag(std::string flag, mmpbsa::MMPBSAState& currState, mmpbsa::MeadInterface& mi);
int parseFlag(std::string flag, mmpbsa::MMPBSAState& currState, mmpbsa::SanderInterface& si);


/**
 * Determines whether or not a filename was provided for the given file type.
 */
bool has_filename(const std::string& filetype, const mmpbsa::MMPBSAState& the_state);

/**
 * Returns a reference to the filename that corresponds to the file type (cf Amber Manual)
 * based on those proved in MMPBSAState. If the file type is not listed in the MMPBSAState
 * object, an exception is thrown. To provent an exception, one could check to see if there
 * is a filename for that file type by using has_filename.
 */
const std::string& get_filename(const std::string& filetype, const mmpbsa::MMPBSAState&) throw (mmpbsa::MMPBSAException);

std::vector<mmpbsa::MMPBSAState> getQueueFile(int argc,char** argv);

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
void checkpoint_out(mmpbsa::MMPBSAState& saveState, mmpbsa_utils::XMLParser& xmlDoc);
void checkpoint_mmpbsa(mmpbsa::MMPBSAState& saveState);
void checkpoint_sander(mmpbsa::MMPBSAState& saveState, mmpbsa::SanderInterface& si);
\
/**
 * Load the last MMPBSA state. Returns true is all of the parameters in the
 * checkpoint file were used.
 *
 */
bool restart_mmpbsa(mmpbsa::MMPBSAState& restartState);
bool restart_sander(mmpbsa::MMPBSAState& restartState, mmpbsa::SanderInterface& si);

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

/**
 * If BOINC and Graphics are being used, communication between the science
 * application and the graphics application is done through the use of shared
 * memory. This method updates the data stored in the shared memory (cf MMPBSA_SHMEM)
 *
 * If BOINC and Graphaics are not used (either or), this function does nothing.
 * 
 */
void update_gshmem();

/**
 * Iterates through the work queue and gives the overall progress of the program
 * as a fraction of one.
 * 
 * @return
 */
double overallFractionDone();

void send_status_message(mmpbsa::SanderInterface& si, double frac_done,
        double checkpoint_cpu_time);

#endif	/* MMPBSA_H */

