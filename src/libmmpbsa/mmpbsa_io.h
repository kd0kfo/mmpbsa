/**
 * Handles IO for MMPBSA classes. 
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
#include "Zipper.h"

#ifdef USE_BOINC
#if defined(_WIN32) || defined(__MINGW_WIN32__)
#include "boinc/boinc_win.h"
#endif
#include "boinc/boinc_api.h"
#endif

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
 * Counts the number of snap shots in the given file.
 * 
 * @param trajFile
 * @param natoms
 * @param isPeriodic
 * @return
 */
size_t count_snapshots(std::fstream& trajFile,const size_t& natoms, bool isPeriodic);

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
std::string getNextLine(std::iostream& file) throw (mmpbsa::MMPBSAException);


/**
 * Reads radii data from a DelPhi file and loads it into the provided maps.
 * If data previously exists in the map that corresponds to a key in the DelPhi
 * file, the data will be overwritten by the DelPhi data.
 * 
 * @param theFile
 * @param radii
 * @param residues
 */
void read_siz_file(std::iostream& theFile,
        std::map<std::string,float>& radii, std::map<std::string,std::string>& residues);


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
template <class T> bool loadValarray(std::iostream& dataFile,
        std::valarray<T>& dataArray, const size_t& arrayLength, const size_t& width,
        const size_t& numberOfColumns);

template <> bool loadValarray<std::string>(std::iostream& dataFile,
            std::valarray<std::string>& dataArray, const size_t& arrayLength, const size_t& width,
            const size_t& numberOfColumns);


/**
 * If BOINC is used, this function will resolve the file name within the BOINC tree
 * and returns the value the boinc function returns.
 * If BOINC is not used, then resolved_name = unresolved_name and zero is returned.
 *
 * @param resolvedFilename
 * @param unresolvedFilename
 */
int resolve_filename(const std::string& unresolvedFilename, std::string& resolvedFilename);

/**
 * If BOINC is used, this function will resolve the file name within the BOINC tree
 * and returns the value the boinc function returns.
 * If BOINC is not used, then resolved_name = unresolved_name and zero is returned.
 *
 * @param resolvedFilename
 * @param unresolvedFilename
 */
int resolve_filename(const char* unresolvedFilename, char* resolvedFilename,  int length);


/**
 * Use stringstream buffers to convert a string representation of a number into the
 * required format. Upon failure a SanderIOException is thrown.
 *
 * @param word string representation of the number.
 * @param data int reference to the variable to be parse.
 * @throws mmpbsa::SanderIOException
 */
void parseNumber(const std::string& word,int& data) throw (mmpbsa::SanderIOException);

/**
 * Use stringstream buffers to convert a string representation of a number into the
 * required format. Upon failure a SanderIOException is thrown.
 *
 * @param word string representation of the number.
 * @param data mmpbsa_t reference to the variable to be parse.
 * @throws mmpbsa::SanderIOException
 */
void parseNumber(const std::string& word, mmpbsa_t& data) throw (mmpbsa::SanderIOException);

/**
 * Use stringstream buffers to convert a string representation of a number into the
 * required format. Upon failure a SanderIOException is thrown.
 *
 * @param word string representation of the number.
 * @param data size_t reference to the variable to be parse.
 * @throws mmpbsa::SanderIOException
 */
void parseNumber(const std::string& word,size_t& data) throw (mmpbsa::SanderIOException);

template <class T> std::ostream& write_snapshot(std::ostream& the_stream,const std::valarray<T>& dataArray,const std::string& ifbox_data);

std::iostream& smart_write(std::iostream& dest, std::iostream& source, const std::string* filename = 0);
 std::iostream& smart_write(std::iostream& dest, const char* source, const size_t& buffer_size, const std::string* filename = 0);

}//end namespace mmpbsa_io


#endif	//SANDERIO_H



