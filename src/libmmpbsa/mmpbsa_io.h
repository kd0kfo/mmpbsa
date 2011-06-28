/**
 * @namespace mmpbsa
 * @brief Input Output functions used by mmpbsa.
 *
 * Created by David Coss <David.Coss@stjude.org> 2010
 */


#ifndef SANDERIO_H
#define	SANDERIO_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <valarray>
#include <iostream>
#include <map>

#include "globals.h"
#include "structs.h"
#include "mmpbsa_exceptions.h"
#include "SanderParm.h"
#include "Vector.h"

#ifdef USE_GROMACS
#include "GromacsReader.h"
#endif

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
std::string get_traj_title(std::istream& trajFile);

/**
 * Counts the number of snap shots in the given file.
 * 
 * @param trajFile
 * @param natoms
 * @param isPeriodic
 * @return
 */
size_t count_snapshots(std::iostream& trajFile,const size_t& natoms, bool isPeriodic);

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
bool get_next_snap(std::iostream& trajFile, std::valarray<mmpbsa::Vector>& snapshot,
    const size_t& natoms,bool isPeriodic = false);

/**
 * DEPRECATED!
 *
 * Skips the next snapshot, but ensures that it had proper data.
 * 
 * @param trajFile
 * @param natoms
 * @param isPeriodic
 */
void skip_next_snap(std::iostream& trajFile, const size_t& natoms,
        bool isPeriodic = false);

/**
 * Reads the next line of a file and returns it as a string.
 * MMPBSAException is thrown if the file cannot be read.
 * 
 * @param file
 * @return 
 */
std::string getNextLine(std::istream& file) throw (mmpbsa::MMPBSAException);

/**
 * Retrieves the next snapshot of the trajectory described by traj.
 *
 * Increments the trajectory_t data field, curr_snap.
 */
bool get_next_snap(mmpbsa_io::trajectory_t& traj, std::valarray<mmpbsa::Vector>& snapshot);

/**
 * Skips ahead (or behind) to the specified snapshot/frame of the
 * trajectory. No data is actually retrieved.
 *
 * snap_pos is one-indexed
 */
void seek(mmpbsa_io::trajectory_t& traj, size_t snap_pos);
void seek(std::iostream& stream, size_t natoms, int ifbox, size_t snap_pos);
/**
 * Initializes the trajectory structure, trajectory_t
 *
 * @see init(mmpbsa_io::trajectory_t*)
 */
void default_trajectory(mmpbsa_io::trajectory_t& traj);

/**
 * Looks into the specified file and determines which type of trajectory
 * it is, e.g. Sander versus Gromacs, and setups a trajectory_t structure
 * for it.
 *
 * Optionally, Sander data can be permanently stored in memory. This is *not*
 * recommended for a large number of snapshots. This is mostly provied for
 * the use compressed data.
 */
mmpbsa_io::trajectory_t open_trajectory(const std::string& filename,const bool& should_remain_in_memory = false);

/**
 * Destructor for a trajectory_t structure. Should be called for each
 * trajectory_t structure, otherwise memory leaks could occur.
 *
 * @see destory(mmpbsa_io::trajectory_t*)
 */
void destroy_trajectory(mmpbsa_io::trajectory_t& traj);

/**
 * Determines whether or not the program has reached the end of
 * the trajectory.
 */
bool eof(trajectory_t& traj);

/**
 * Retrieves the title of the trajectory.
 */
std::string get_traj_title(mmpbsa_io::trajectory_t& traj);


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
template <> bool loadValarray<mmpbsa::Vector>(std::iostream& dataFile,
            std::valarray<mmpbsa::Vector>& dataArray, const size_t& arrayLength, const size_t& width,
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

/**
 * Based on the filename, smart_write decides whether or not compression should
 * be used and if so calls the appropriate function. If neither of the streams
 * represent a file, set the filename pointer to 0 (null).
 * 
 * iostream is used to extend the possible object to either file or string streams.
 *
 * @param std::iostream& dest Destination stream. 
 * @param std::iostream& source source stream. 
 * @param const std::string* pointer to filename string (default = NULL)
 * @return std::iostream& reference to destination stream.
 */
std::iostream& smart_write(std::iostream& dest, std::iostream& source, const std::string* filename = 0);

/**
 * Based on the filename, smart_write decides whether or not compression should
 * be used and if so calls the appropriate function. Uses the same procedure
 * as the stream-source smart_write function, except a character array is 
 * written. 
 * If neither source nor destination represent a file, set the filename pointer to 0 (null).
 * 
 * iostream is used to extend the possible object to either file or string streams.
 *
 * @param std::iostream& dest Destination stream. 
 * @param const char* source buffer (need not be null terminated)
 * @param const size_t& buffer_size Size of the buffer, in units of sizeof(char)
 * @param const std::string* pointer to filename string (default = NULL)
 * @return std::iostream& reference to destination stream.
 */
std::iostream& smart_write(std::iostream& dest, const char* source, const size_t& buffer_size, const std::string* filename = 0);

/**
 * Based on the filename, smart_read decides whether or not compression should
 * be used and if so calls the appropriate function. Uses the same procedure
 * as the stream-destination smart_read function, except data is placed in  
 * a char array.
 * If neither source nor destination represent a file, set the filename pointer to 0 (null).
 * 
 * iostream is used to extend the possible object to either file or string streams.
 *
 * @param const char* destination buffer (need not be null terminated, i.e. binary data)
 * @param std::iostream& dest source stream. 
 * @param const std::string* pointer to filename string (default = NULL)
 * @return size_t size of destination buffer.
 */
size_t smart_read(char** dest, std::iostream& source, const std::string* filename = 0);

/**
 * Based on the filename, smart_read decides whether or not compression should
 * be used and if so calls the appropriate function. 
 *
 * If neither source nor destination represent a file, set the filename pointer to 0 (null).
 * 
 * iostream is used to extend the possible object to either file or string streams.
 *
 * @param std::iostream& destination stream.
 * @param std::iostream& dest source stream. 
 * @param const std::string* pointer to filename string (default = NULL)
 * @return size_t size of destination buffer.
 */
std::iostream& smart_read(std::iostream& dest, std::iostream& source, const std::string* filename = 0);

std::string pdbPad(const int& neededDigits,const int& currentNumber);

}//end namespace mmpbsa_io

/**
 * Initializes a mmpbsa_io::trajectory_t structure
 */
void init(mmpbsa_io::trajectory_t* traj);

/**
 * Initializes a mmpbsa_io::trajectory_t structure.
 * Needs to be called for each trajectory_t structure.
 * Otherwise, memory could be leaked.
 */
void destroy(mmpbsa_io::trajectory_t* traj);

/**
 * Creates a PDB for the lsit of atoms in the given forcefield and
 * streams it to the output stream.
 */
std::ostream& streamPDB(std::ostream& theStream, const std::vector<mmpbsa::atom_t>& atoms,const mmpbsa::forcefield_t& ff, const std::valarray<mmpbsa::Vector>& crds) throw (mmpbsa::MMPBSAException);

#endif	//SANDERIO_H



