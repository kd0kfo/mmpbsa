/**
 * Zipper
 *
 * Compress routines. Can use a combination of tar and zlib archiving and compressions.
 *
 * zip and unzip compression routines adapted from public source code from zlib's website, www.zlib.net
 * These rouintes produce *gzip* files.
 *
 * Created by David Coss, 2007
 *
 */

#ifndef ZIPPER_H_
#define ZIPPER_H_

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <zlib.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mmpbsa_exceptions.h"

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#define SET_BINARY_MODE(file)
#endif

#define CHUNK 16384
#define ZIPPER_WINDOW_BITS 31
#define TAR_BLOCK_SIZE 512

static const int DEFAULT_COMPRESSION = Z_DEFAULT_COMPRESSION;

namespace mmpbsa_utils{

class Zipper {
public:
	/**
	 * Compress from file source to file dest until EOF on source.
	 * def() returns Z_OK on success, Z_MEM_ERROR if memory could not be
	 * allocated for processing, Z_STREAM_ERROR if an invalid compression
	 * level is supplied, Z_VERSION_ERROR if the version of zlib.h and the
	 * version of the library linked do not match, or Z_ERRNO if there is
	 * an error reading or writing the files.
	 */
	static int fzip(FILE *source, FILE *dest, int level);
	static int zip(const std::string& source_filename, const std::string& dest_filename);

	/**
	 * Decompress from file source to file dest until stream ends or EOF.
	 * inf() returns Z_OK on success, Z_MEM_ERROR if memory could not be
	 * allocated for processing, Z_DATA_ERROR if the deflate data is
	 * invalid or incomplete, Z_VERSION_ERROR if the version of zlib.h and
	 * the version of the library linked do not match, or Z_ERRNO if there
	 * is an error reading or writing the files.
	 */
	static int funzip(FILE *source, FILE *dest);
	static int unzip(const std::string& source_filename, const std::string& dest_filename);

	static void zerr(int ret);

	static void pad_tarfile(const char* data_buffer, FILE* out_file);
	static void create_header(char* header, const char* data, const size_t& data_size, const struct stat& file_stat, const char* data_filename);
	static struct stat default_stat(const size_t& file_size);

	/**
	 * Takes a list of files and creates a tar file.
	 * The directory of the files is specified by dir. The prefix of the files,
	 * i.e. the directory of the files after the tar file is unpacked, is
	 * specified by prefix.
	 *
	 * NOTE: Does not close out_file. Only writes to it.
	 */
	static void tar(const std::vector<std::string>& filenames,FILE* out_file,
			const std::string& dir,const std::string& prefix);

	/**
	 * Returns a pointer to a stringstream containing the data corresponding to filename.
	 *
	 * If filename is not found in the tar headers, a null pointer is returned.
	 */
	static std::stringstream* funtar(FILE* in_file,const std::string& filename);

	/**
	 * Internal routine that writes a number to a char* buffer in octal radix.
	 * padding (to the left) with zeros so that the width is equal to the
	 * provided width.
	 */
	static void write_oct(char* buffer, const int& number, const int& width);
	static void write_dec(char* buffer, const int& number, const int& width);

	/**
	 * tar checksum.
	 * Sum of all char values in the header as unsigned bytes, when the
	 * checksum position is filled with space characters (0x20).
	 */
	static int checksum(const char* header,const size_t& header_size);


};

}//namespace mmpbsa_utils

namespace mmpbsa{

class ZipperException : public mmpbsa::MMPBSAException
{
public:
	ZipperException(const std::string& error) : mmpbsa::MMPBSAException( error){}

	ZipperException(const std::string& error, const mmpbsa::MMPBSAErrorTypes& errorType)
	: mmpbsa::MMPBSAException(error,errorType){}
	ZipperException(const std::ostringstream& error, const mmpbsa::MMPBSAErrorTypes& errorType)
	: mmpbsa::MMPBSAException(error,errorType){}
	ZipperException(const std::ostringstream& error)
		: mmpbsa::MMPBSAException(error,UNKNOWN_ERROR){}

	const char* identifier(){return "Compression Error";}

};

}//namespace mmpbsa

#endif /* ZIPPER_H_ */
