#ifndef MMPBSA_H
#define	MMPBSA_H

#include <cstdlib>
#include <iostream>
#include <valarray>
#include <fstream>

#include "mmpbsa_exceptions.h"
#include "mmpbsa_io.h"
#include "EnergyInfo.h"
#include "Energy.h"
#include "MeadInterface.h"

#include "MEAD/FinDiffMethod.h"

/**
 * Pulls everything together to perform the MMPBSA calculations.
 * Calls the parseArgs to load information from the command line.
 * 
 * @param argc
 * @param argv
 * @return
 */
int realDeal(int argc, char** argv);
void printSnapshot(const EMap& complexEMap, const EMap& receptorEMap,
        const EMap& ligandEMap, std::fstream& outFile);

/**
 * Takes the command line arguments and performs the necessary functions.
 * 
 * @param argc
 * @param argv
 * @param mi
 */
void parseArgs(int argc, char** argv, MeadInterface& mi);

/**
 * When a command line argument provides data or information(to the right of an
 * equal sign), this method decides what to do with it.
 *
 * @param arg
 * @param mi
 */
void parseParameter(std::string arg, MeadInterface& mi);

/**
 * When a command line argument toggles a program flag, this method makes the
 * proper changes to the program.
 * 
 * @param flag
 * @param mi
 */
void parseFlag(std::string flag, MeadInterface& mi);

/**
 * When a parameter has a list of variables (separated by commas), this method
 * tokenizes that list and loades it into the array, in left to right order.
 * 
 * @param values
 * @param array
 */
void loadListArg(const std::string& values,std::vector<size_t>& array);

/**
 * Creates a string about the usage of the program, listing the parameters and
 * flags. This string should be sent to STDOUT, if "--help" is provided as a flag.
 * 
 * @return
 */
std::string helpString();



class MMPBSAState
{
public:
    std::vector<size_t> receptorStartPos;///<Starting positions for receptors. End positions are deduced from the parmtop file.
    std::vector<size_t> ligandStartPos;///<Starting positions for ligands. End positions are deduced from the parmtop file.
    std::vector<size_t> snapList;///<List of snapshots to be used in the calculations. If the list is empty, all snapshots are used.

    std::fstream prmtopFile;
    std::fstream trajFile;
    std::fstream radiiFile;
    std::fstream outputFile;

    std::string checkpointFilename;

    size_t currentSnap;
    enum MOLECULE{COMPLEX,RECEPTOR,LIGAND,END_OF_MOLECULES};///<List parts of the complex on which to perform MMPBSA. END_OF_MOLECUES signifies the snapshot is finished.
    MOLECULE currentMolecule;

    bool trustPrmtop;///<Flag to indicate if the sanity check of the SanderParm object should be ignored. This is not suggested, but if the sanity check fails and one *does* believe it should work, this is provided as a work around, for whatever reason might arise.
    /**
     * Stores variables needed to restart the program. This is needed for running
     * on systems where the program may be closed and restarted again, i.e. BOINC Grid
     */
    MMPBSAState();
    MMPBSAState(const MMPBSAState& orig);
    virtual ~MMPBSAState(){}

    /**
     * Saves the state of the MMPBSA Program. If BOINC is being used, BOINC checkpointing
     * is performed.
     */
    void checkpoint_out();

    /**
     * Load the last MMPBSA state. Returns true is all of the parameters in the
     * checkpoint file were used.
     * 
     */
    bool checkpoint_in(const std::string& fileName);
    bool checkpoint_in(){return checkpoint_in(this->checkpointFilename);}

    /**
     * Flushes all of the buffers within the MMPBSA state (if they're open at all;
     * MMPBSAState tries not to be stingy with files).
     */
    void flush();

    /**
     * Closes all of the buffers within the MMPBSA state.
     */
    void close();

    MMPBSAState& operator=(const MMPBSAState& orig);
private:
    size_t checkpointCounter;
};

MMPBSAState currState;

#endif	/* MMPBSA_H */

