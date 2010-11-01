
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
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
    filename_map = orig.filename_map;
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
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
    filename_map = orig.filename_map;

    return *this;
}
