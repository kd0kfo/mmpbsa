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

int realDeal(int argc, char** argv);
void printSnapshot(const EMap& complexEMap, const EMap& receptorEMap,
        const EMap& ligandEMap, std::fstream& outFile);
int testsubmain(int argc, char** argv);
int testenergy(int argc, char** argv);
int testPolyAT(int argc, char** argv);
void parseArgs(int argc, char** argv, MeadInterface& mi);
void parseParameter(std::string arg, MeadInterface& mi);
void parseFlag(std::string flag, MeadInterface& mi);
void loadListArg(const std::string& values,std::vector<size_t>& array);
std::string helpString();

std::vector<size_t> receptorStartPos;
std::vector<size_t> ligandStartPos;
std::vector<size_t> snapList;


std::fstream myOutput;
std::fstream prmtopFile;
std::fstream trajFile;
std::fstream radiiFile;
std::fstream outputFile;

#endif	/* MMPBSA_H */

