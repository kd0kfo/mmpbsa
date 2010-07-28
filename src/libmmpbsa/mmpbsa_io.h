/**
 * SanderIO
 *
 * Handles IO for sander parameter files. Parameters are stored in the
 *      "SanderParm" class.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */


#ifndef SANDERIO_H
#define	SANDERIO_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//Standard Includes
#include <string>
#include <cstring>
#include <fstream>
#include <valarray>
#include <streambuf>
#include <sstream>
#include <iostream>

//project specific stuff
#include "mmpbsa_utils.h"
#include "mmpbsa_exceptions.h"
#include "SanderParm.h"

namespace mmpbsa_io{

/**
 * Reads the Coordinates or Velocities from the given file, using the first
 * natoms*3 values. The title (first line value) is returned as a std::string
 *
 * @param crdFile
 * @param crds
 * @return title
 */
std::string read_crds(std::fstream& crdFile, std::valarray<mmpbsa_t>& crds);

/**
 * Writes the provided coordinate data to a file. Each line contains 6 columns
 * with 12.7f formatting.
 *
 * @param fileName
 * @param crds
 * @param title
 */
void write_crds(const char* fileName,const std::valarray<mmpbsa_t>& crds,
    const char* title = "");

/**
 * Moves the current read line of the file stream to the beginning, increments
 * to the first snapshot and returns the title.
 * 
 * @param trajFile
 * @return title
 */
std::string get_traj_title(std::fstream& trajFile);

/**
 * Gets the next snapshot from the provided trajectory file. The snapshot data
 * is loaded into the provided snapshot valarray, overwriting data and resizing,
 * if necessary. True is returned in the snapshot is read. False otherwise, ie
 * trajFile is at EOF.
 * 
 * @param trajFile
 * @param snapshot
 * @param natoms
 * @return
 */
bool get_next_snap(std::fstream& trajFile, std::valarray<mmpbsa_t>& snapshot,
    const size_t& natoms,bool isPeriodic = false);

/**
 * Skips the next snapshot, but ensures that it had proper data.
 * 
 * @param trajFile
 * @param natoms
 * @param isPeriodic
 */
void skip_next_snap(std::fstream& trajFile, const size_t& natoms,
        bool isPeriodic = false);

/**
 * Reads the next line of a file and returns it as a string.
 * MMPBSAException is thrown if the file cannot be read.
 * 
 * @param file
 * @return 
 */
std::string getNextLine(std::fstream& file) throw (mmpbsa::MMPBSAException);


/**
 * Reads radii data from a DelPhi file and loads it into the provided maps.
 * If data previously exists in the map that corresponds to a key in the DelPhi
 * file, the data will be overwritten by the DelPhi data.
 * 
 * @param theFile
 * @param radii
 * @param residues
 */
void read_siz_file(std::fstream& theFile,
        std::map<std::string,mmpbsa_t>& radii, std::map<std::string,std::string>& residues);

/**
 * Opens a file using the provided fstream.
 *
 * If the program is compiled with the BOINC API, the return value of boinc_resolve_filename
 * is returned. Otherwise, the return values are 0 for success and 1 for failure.
 *
 * @param dataFile
 * @param dataArray
 * @param arrayLength
 * @param width
 * @param numberOfColumns
 * @return
 */
int fileopen(const char* filename, const std::ios::openmode& mode,
        std::fstream& file);


/**
 * Takes a given file and loads data into the given valarray, overwriting data
 * if it already exists in valarray. The size of the valarray is set to arrayLength,
 * if it is not already. The width of each column (including whitespace) is 
 * equal to width. True is returned in the snapshot is read. False otherwise, ie
 * trajFile is at EOF.
 * 
 * @param dataFile
 * @param dataArray
 * @param arrayLength
 * @param width
 * @return 
 */
template <class T> bool loadValarray(std::fstream& dataFile,
        std::valarray<T>& dataArray, const size_t& arrayLength, const size_t& width,
        const size_t& numberOfColumns);

template <> bool loadValarray<std::string>(std::fstream& dataFile,
            std::valarray<std::string>& dataArray, const size_t& arrayLength, const size_t& width,
            const size_t& numberOfColumns);

/**
     * Gets the next line with data, ie empty, whitespace lines are ignored.
     * @param file
     * @return
     */

}//end namespace mmpbsa_io



#endif	//SANDERIO_H



