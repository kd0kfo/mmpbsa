#include "mmpbsa_io.h"

std::string mmpbsa_io::read_crds(std::fstream& crdFile, std::valarray<mmpbsa_t>& crds)
{
    using std::string;
    using namespace mmpbsa_utils;
    
    if(!crdFile.good())
        throw mmpbsa::SanderIOException("Cannot open coordinate file",mmpbsa::FILE_IO_ERROR);

    string title = getNextLine(crdFile);
    string strNatoms = getNextLine(crdFile);
    strNatoms =trimString(strNatoms);
    size_t natoms = 0;
    std::istringstream buff(strNatoms);
    buff >> natoms;

    if(!loadValarray(crdFile,crds,natoms*3,12,8))
        throw mmpbsa::SanderIOException("Coordinate file is too short.",mmpbsa::FILE_IO_ERROR);

    return title;
}

void mmpbsa_io::write_crds(const char* fileName,const std::valarray<mmpbsa_t>& crds,
    const char* title)
{
    using std::valarray;
    using std::slice;

    if(crds.size() % 3 != 0)
        throw mmpbsa::SanderIOException("The number of elements in the coordinate array "
                "must be a multiple of 3, ie 3-dimensions.",mmpbsa::DATA_FORMAT_ERROR);

    size_t natoms = size_t(crds.size()/3);

    std::fstream outFile(fileName,std::ios::out);
    
    if(!outFile.good())
    {
        std::ostringstream error;
        error << "Could not open: " << fileName;
        throw mmpbsa::SanderIOException(error,mmpbsa::FILE_IO_ERROR);
    }

    outFile << title << std::endl;
    outFile << std::setprecision(5) << natoms << std::endl;

    size_t m;
    //save data in rows of 6
    outFile << std::setprecision(7);
    outFile.width(8);
    for(m = 0;m<crds.size() - 6;m+=6)
    {
        valarray<mmpbsa_t> row = crds[slice(m,6,1)];//m-th row

        for(size_t i = 0;i<5;i++)
        {
            outFile << row[i] << " ";
        }
        outFile << row[5] << std::endl;//no need for " " after the last entry
    }

    //save the possibly incomplete last row.
    if(m<crds.size())
    {
        for(;m<crds.size()-1;m++)
        {
            outFile << crds[m] << " ";
        }
        outFile << crds[crds.size()-1] << std::endl;
    }

    outFile.close();
}


std::string mmpbsa_io::get_traj_title(std::iostream& trajFile)
{
    trajFile.seekg(0,std::ios::beg);
    return getNextLine(trajFile);
}

std::string mmpbsa_io::getNextLine(std::iostream& file) throw (mmpbsa::MMPBSAException)
{
    if(!file.good())
        throw mmpbsa::MMPBSAException("Could not read from file");

    std::string returnMe;
    getline(file,returnMe);
    return returnMe;
}



bool mmpbsa_io::get_next_snap(std::iostream& trajFile, std::valarray<mmpbsa_t>& snapshot,
    const size_t& natoms,bool isPeriodic)
{
    bool returnMe = loadValarray(trajFile,snapshot,natoms*3,8,10);
    if(isPeriodic)
        getNextLine(trajFile);//ignoring periodic box information
    return returnMe;
}

size_t mmpbsa_io::count_snapshots(std::iostream& trajFile,const size_t& natoms, bool isPeriodic)
{
    get_traj_title(trajFile);
    size_t snapcount = 0;
    try
    {
        while(!trajFile.eof())
        {
            skip_next_snap(trajFile,natoms,isPeriodic);
            snapcount++;
        }
    }
    catch(mmpbsa::SanderIOException sioe)
    {
        if(sioe.getErrType() == mmpbsa::UNEXPECTED_EOF)
            return snapcount;
        throw sioe;
    }
    return snapcount;
}

void mmpbsa_io::skip_next_snap(std::iostream& trajFile, const size_t& natoms, bool isPeriodic)
{
    throw mmpbsa::MMPBSAException("mmpbsa_io::skip_next_snap: DEPRECATED!");
}

void mmpbsa_io::parseNumber(const std::string& word,int& intData) throw (mmpbsa::SanderIOException)
{
	std::istringstream dataBuffer(word);
	dataBuffer >> intData;
	if(dataBuffer.fail())
	{
		throw mmpbsa::SanderIOException("parseNumber expected an integer but received \""
				+ word + "\"",mmpbsa::DATA_FORMAT_ERROR);
	}
}

void mmpbsa_io::parseNumber(const std::string& word,size_t& intData) throw (mmpbsa::SanderIOException)
{
	std::istringstream dataBuffer(word);
	dataBuffer >> intData;
	if(dataBuffer.fail())
	{
		throw mmpbsa::SanderIOException("parseNumber expected an integer but received \""
				+ word + "\"",mmpbsa::DATA_FORMAT_ERROR);
	}
}

void mmpbsa_io::parseNumber(const std::string& word,mmpbsa_t& returnData) throw (mmpbsa::SanderIOException)
{
	std::istringstream dataBuffer(word);
	dataBuffer >> returnData;
	if(dataBuffer.fail())
	{
		throw mmpbsa::SanderIOException("parseNumber expected an mmpbsa data type but received \""
				+ word + "\"",mmpbsa::DATA_FORMAT_ERROR);
	}
}


template <class T> bool mmpbsa_io::loadValarray(std::iostream& dataFile,
        std::valarray<T>& dataArray, const size_t& arrayLength, const size_t& width,
        const size_t& numberOfColumns)
{
    using std::string;

    //If the length is zero, there is no data, which will correspond to a blank
    //line in the parmtop file. Pop that line and return (true);
    if(arrayLength == 0)
    {
        getNextLine(dataFile);
        return true;
    }

    if(dataFile.eof())
        return false;

    if(dataArray.size() != arrayLength)
        dataArray.resize(arrayLength);

    size_t lineIndex = 0;
    size_t dataIndex = 0;

    for(;dataIndex<arrayLength;)
    {
        if(dataFile.eof())
            throw mmpbsa::SanderIOException("Data file ended in the middle of the "
                    "data.",mmpbsa::UNEXPECTED_EOF);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
        	std::cerr << "Data file contains a short line. "
        			"Lines must be at least 36 characters, but line #"
        			<< lineIndex+1 << " is only " << currentLine.size()
        			<< " characters long." << std::endl;
        }

        //tokenize line into data. put data into valarray.
        while(currentLine.size() >= width)
        {
        	parseNumber(currentLine.substr(0,width),dataArray[dataIndex++]);
            currentLine.erase(0,width);
        }

        lineIndex++;
    }

    return true;

}

template <> bool mmpbsa_io::loadValarray<std::string>(std::iostream& dataFile,
        std::valarray<std::string>& dataArray, const size_t& arrayLength, const size_t& width,
        const size_t& numberOfColumns)
{
    using std::string;

    if(dataFile.eof())
        return false;

    //If the length is zero, there is no data, which will correspond to a blank
    //line in the parmtop file. Pop that line and return (true);
    if(arrayLength == 0)
    {
        getNextLine(dataFile);
        return true;
    }

    if(dataArray.size() != arrayLength)
        dataArray.resize(arrayLength,"");

    size_t lineIndex = 0;
    size_t dataIndex = 0;

    for(;dataIndex<arrayLength;)
    {
        if(dataFile.eof())
            throw mmpbsa::SanderIOException("Data file ended in the middle of the "
                    "data.",mmpbsa::BROKEN_TRAJECTORY_FILE);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
            std::cerr << "Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #"
                    << lineIndex+1 << " is only "
                    << currentLine.size() <<  " characters long." << std::endl;
        }

        //tokenize line into data. put data into valarray.
        while(currentLine.size() > 0)
        {
            dataArray[dataIndex++] = currentLine.substr(0,width);
            currentLine.erase(0,width);
        }

        lineIndex++;
    }


    return true;
}

template <class T> std::ostream& mmpbsa_io::write_snapshot(std::ostream& the_stream,const std::valarray<T>& dataArray,const std::string& ifbox_data)
{
	//If the length is zero, there is no data, which will correspond to a blank
	//line in the parmtop file. Pop that line and return (true);
	if(dataArray == 0)
	{
		std::cerr << "Warning: write_snapshot was called without data." << std::endl;
		return the_stream;
	}

	size_t arrayLength = dataArray.size();
	if(arrayLength % 3 != 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::write_snapshot: Given coordinate (trajectory) data that is not 3-dimensional.",mmpbsa::DATA_FORMAT_ERROR);

	for(size_t dataIndex = 0;dataIndex<arrayLength;dataIndex++)
	{
		if(dataIndex % 10 == 0)
		{
			the_stream << std::endl;
		}

		the_stream << std::setw(8) << dataArray[dataIndex];
	}
	if(ifbox_data.size() != 0)
		the_stream << std::endl << ifbox_data;
	the_stream << std::endl;
	return the_stream;
}

void mmpbsa_io::read_siz_file(std::iostream& theFile,
        std::map<std::string,float>& radii, std::map<std::string,std::string>& residues)
{
    using mmpbsa_utils::trimString;
    using mmpbsa_utils::toUpperCase;
    
    if(!theFile.good())
        throw mmpbsa::MMPBSAException("Could not open SIZ file.",mmpbsa::FILE_IO_ERROR);

    std::string currLine;
    std::string atomName;
    std::string residue;
    std::string data;
    mmpbsa_t fData;
    size_t lineNumber = 0;
    while(theFile.good())
    {
        currLine = getNextLine(theFile);
        if(trimString(currLine) == "")//there may be a blank line at the end. Ignore blank lines
            continue;
        lineNumber++;
        if(currLine[0] == '!')//comments begin with "!"
            continue;
        if(currLine.substr(0,5) == "atom_")//the format line in the file will begin with "atom_". Perhaps later dynamically read format??
            continue;
        if(currLine.size() < 9)//data begins after atomname(6chars) and residue name (3chars)
        {
            std::ostringstream error;
            error << "Improperly formatted SIZ file: Short line at " << lineNumber;
            throw mmpbsa::MMPBSAException(error,mmpbsa::FILE_IO_ERROR);
        }

        atomName = toUpperCase(trimString(currLine.substr(0,6)));
        switch(*(atomName.begin()))
        {
            //Atom numbering in DelPhi and Parmtop files are reversed.
            case '1': case '2': case '3': case '4': case '5':
                atomName = atomName.substr(1).append(atomName.substr(0,1));
                break;
            default:
                break;
        }
        residue = toUpperCase(trimString(currLine.substr(6,3)));
        std::istringstream dataStream(trimString(currLine.substr(9)));
        dataStream >> MMPBSA_FORMAT >> fData;

        radii[atomName] = fData;
        residues[atomName] = residue;
    }

}


std::iostream& mmpbsa_io::smart_write(std::iostream& dest, std::iostream& source, const std::string* filename)
{

	if(!dest.good())
			throw mmpbsa::MMPBSAException("smart_write: cannot write to stream.");

	if(filename == 0)
	{
		dest << source.rdbuf();
		return dest;
	}

	bool should_gzip = (filename->find(".gz") != std::string::npos || filename->find(".tgz") != std::string::npos);
	bool should_tar = (filename->find(".tgz") != std::string::npos || filename->find(".tar") != std::string::npos);
	if(!should_gzip && !should_tar)//for efficiency's sake, if there is no compression, don't use buffers.
	{
		dest << source.rdbuf();
		return dest;
	}


	size_t buffer_size;
	char *buffer;
	source.seekg(0,std::ios::end);
	buffer_size = source.tellg();
	source.seekg(0,std::ios::beg);
	buffer = new char[buffer_size];

	smart_write(dest,buffer,buffer_size,filename);

	delete [] buffer;
	return dest;
}

std::iostream& mmpbsa_io::smart_write(std::iostream& dest, const char* source,const size_t& buffer_size, const std::string* filename)
{
	if(!dest.good())
		throw mmpbsa::MMPBSAException("smart_write: cannot write to stream.");

	if(source == 0)
		throw mmpbsa::MMPBSAException("smart_write: source buffer is a null pointer.",mmpbsa::NULL_POINTER);

#ifdef USE_GZIP
	using mmpbsa_utils::Zipper;
	FILE *tempfile,*out_file = tmpfile();

	bool should_gzip = (filename->find(".gz") != std::string::npos || filename->find(".tgz") != std::string::npos);
	bool should_tar = (filename->find(".tgz") != std::string::npos || filename->find(".tar") != std::string::npos);
	if(should_gzip || should_tar)
	{
		char header[TAR_BLOCK_SIZE],*buffer;
		std::string decomp_filename = *filename;
		struct stat the_stat = Zipper::default_stat(buffer_size);
		size_t gzip_size;
		tempfile = tmpfile();

		//Tar file
		if(should_tar)
		{
			if(decomp_filename.find(".tgz") != std::string::npos)
				decomp_filename = decomp_filename.substr(0,decomp_filename.find(".tgz"));
			else if(decomp_filename.find(".tar") != std::string::npos)
				decomp_filename = decomp_filename.substr(0,decomp_filename.find(".tar"));
			Zipper::create_header(header,source,the_stat.st_size,the_stat,decomp_filename.c_str());
			fwrite(header,1,TAR_BLOCK_SIZE,tempfile);
			fwrite(source,1,the_stat.st_size,tempfile);
			Zipper::pad_tarfile(source,tempfile);
			rewind(tempfile);
		}

		//gzip file (or intermediate tar file)
		if(should_gzip)
		{
			Zipper::fzip(tempfile,out_file,Z_BEST_SPEED);
			fclose(tempfile);
		}
		else
		{
			fclose(out_file);
			out_file = tempfile;
		}

		fseek(out_file,0,SEEK_END);
		gzip_size = ftell(out_file);
		buffer = new char[gzip_size];
		fseek(out_file,0,SEEK_SET);
		fread(buffer,sizeof(char),gzip_size,out_file);
		dest.write(buffer,gzip_size);
		delete [] buffer;
	}
	else//in this case, no tar or gzip
		dest << source;

	fclose(out_file);
#else
	dest << source;
#endif

	return dest;
}

size_t mmpbsa_io::smart_read(char** dest, std::iostream& source, const std::string* filename)
{
	size_t buffer_size;
#ifdef USE_GZIP
	using mmpbsa_utils::Zipper;
	FILE *tempfile,*in_file = tmpfile();

	char* buffer;
	size_t source_size;

	//push the buffer into the FILE stream to be used by tar and/or gzip
	source.seekg(0,std::ios::end);
	source_size = source.tellg();
	source.seekg(0,std::ios::beg);
	buffer = new char[source_size];
	source.read(buffer,source_size);
	fwrite(buffer,sizeof(char),source_size,in_file);
	rewind(in_file);

	//Without an extension, do not try to decompress.
	if(filename == 0)
	{
		fclose(in_file);
		*dest = buffer;
		return source_size;
	}

	bool should_gzip = (filename->find(".gz") != std::string::npos || filename->find(".tgz") != std::string::npos);
	bool should_tar = (filename->find(".tgz") != std::string::npos || filename->find(".tar") != std::string::npos);
	if(should_gzip || should_tar)
	{
		//gzip file (or intermediate tar file)
		if(should_gzip)
		{
			tempfile = tmpfile();
			int unzip_result = Zipper::funzip(in_file,tempfile);
			if(unzip_result != Z_OK)
			{
				std::ostringstream error;
				Zipper::zerr(unzip_result);
				error << "write_mmpbsa_data: Error decompressing " << *filename << " zlib error number " << unzip_result;
				throw mmpbsa::MMPBSAException(error);
			}
			rewind(tempfile);
			fclose(in_file);
		}
		else
			tempfile = in_file;

		//Tar file
		if(should_tar)
		{
			std::string decomp_filename;
			if(filename->find(".tgz") != std::string::npos)
				decomp_filename = filename->substr(0,filename->find(".tgz"));
			else if(filename->find(".tar") != std::string::npos)
				decomp_filename = filename->substr(0,filename->find(".tar"));
			std::stringstream* tarstream = Zipper::funtar(tempfile,decomp_filename);
			if(tarstream == 0)
				throw mmpbsa::ZipperException("mmpbsa_io::smart_read: bad tar data in " + *filename,mmpbsa::FILE_IO_ERROR);
			buffer_size = tarstream->str().size();
			*dest = new char[buffer_size];
			memcpy(*dest,tarstream->str().c_str(),buffer_size);
			delete tarstream;
		}
		else
		{
			fseek(tempfile,0,SEEK_END);
			buffer_size = ftell(tempfile);
			fseek(tempfile,0,SEEK_SET);
			*dest = new char[buffer_size];
			fread(*dest,sizeof(char),buffer_size,tempfile);
		}

		fclose(tempfile);
		delete [] buffer;
		return buffer_size;
	}//end should tar or should gzip

	//If this point is reached, no tar or gzip; just use XMLParser.
	*dest = buffer;
	fclose(in_file);
	return source_size;

#else
	source.seekg(0,std::ios::end);
	buffer_size = source.tellg();
	source.seekg(0,std::ios::beg);
	*dest = new char[buffer_size];
	source.read(*dest,buffer_size);
	return buffer_size;
#endif
}

std::iostream& mmpbsa_io::smart_read(std::iostream& dest, std::iostream& source, const std::string* filename)
{
	if(!dest.good())
		throw mmpbsa::MMPBSAException("smart_read: cannot write to stream.");

	if(!source.good())
			throw mmpbsa::MMPBSAException("smart_read: cannot read input stream.");

	if(filename == 0)
	{
		dest << source.rdbuf();
		return dest;
	}

	bool should_gzip = (filename->find(".gz") != std::string::npos || filename->find(".tgz") != std::string::npos);
	bool should_tar = (filename->find(".tgz") != std::string::npos || filename->find(".tar") != std::string::npos);
	if(!should_gzip && !should_tar)//save time and not use buffers.
	{
		dest << source.rdbuf();
		return dest;
	}

	char *buffer = 0;
	size_t buffer_size = smart_read(&buffer,source,filename);
	if(buffer_size != 0 && buffer != 0)
	  dest.write(buffer,buffer_size);
	delete [] buffer;
	return dest;
}

int mmpbsa_io::resolve_filename(const std::string& unresolvedFilename, std::string& resolvedFilename)
{
#ifdef USE_BOINC
      return boinc_resolve_filename_s(unresolvedFilename.c_str(),resolvedFilename);
#else
    resolvedFilename = unresolvedFilename;
    return 0;
#endif
}

#if !defined(HAVE_STRLCPY) && !defined(USE_BOINC)
//These are needed for start() below. However, if BOINC is not used, it must
//be provided here.
//Copied from str_util.cpp under the terms of the GNU Lesser General Public License.
//See http://boinc.berkeley.edu or http://www.gnu.org for details.
size_t strlcpy(char *dst, const char *src, size_t size) {
    size_t ret = strlen(src);

    if (size) {
        size_t len = (ret >= size) ? size-1 : ret;
        memcpy(dst, src, len);
        dst[len] = '\0';
    }

    return ret;
}
#endif

int mmpbsa_io::resolve_filename(const char* unresolvedFilename, char* resolvedFilename, int length)
{
#ifdef USE_BOINC
    return boinc_resolve_filename(unresolvedFilename, resolvedFilename,
            length);
#else
    if(!unresolvedFilename)
        return 1;
    strlcpy(resolvedFilename,unresolvedFilename,length);
    return 0;
#endif
}


std::string mmpbsa_io::pdbPad(const int& neededDigits,const int& currentNumber)
{
	std::string returnMe = "";
	for(size_t serialPad = 0;serialPad < neededDigits - floor(log(currentNumber)/log(10.0));serialPad++)
					returnMe += " ";
	return returnMe;
}

std::ostream& streamPDB(std::ostream& theStream, const std::vector<mmpbsa::atom_t>& atoms,const mmpbsa::forcefield_t& ff, const std::valarray<mmpbsa_t>& crds) throw (mmpbsa::MMPBSAException)
{
	using mmpbsa_io::pdbPad;
	using namespace mmpbsa;

	if(!theStream.good())
		throw MMPBSAException("Could not write to the stream provided to streamPDB",FILE_IO_ERROR);

	//iterate through list. output using PDB format.
	size_t currRes = 1,currAtom = 0;;
	size_t serialNumber = 1;
	std::vector<atom_t>::const_iterator atom = atoms.begin();
	for(;atom != atoms.end();atom++,currAtom++)
	{
		//record name (cols 1-6)
		theStream << "ATOM  ";

		//serial number (cols 7-11)
		theStream << pdbPad(4,serialNumber);
		theStream << serialNumber++;

		//empty space in column 12
		theStream << " ";

		//atom name (from parm top file) (cols 13-16)
		theStream << atom->name;

		//alt location indicator (col 17)
		theStream << " ";

		//Residue name(cols 18-20)
		std::string resName = "DRC";
		theStream << resName.substr(0,3);

		//empty space column 21 chain id column 22
		theStream << "  ";

		//residue sequence (cols 23 - 26)
		theStream << pdbPad(4,currRes);
		theStream << currRes;

		//icode (col 27) with columns 28-30 empty
		theStream << "    ";

		//Coordinates (cols 31 - 54)
		std::ostringstream coordinateBuffer;
		coordinateBuffer.setf(std::ios::fixed,std::ios_base::floatfield);
		coordinateBuffer.precision(3);
		for(size_t coordIndex = 0;coordIndex < 3;coordIndex++)
		{
			coordinateBuffer.width(8);
			coordinateBuffer << crds[3*currAtom + coordIndex];
			theStream << coordinateBuffer.str();
			coordinateBuffer.str("");
		}

		//occupancy & temp factor
		theStream << "  1.00" << "  0.00";

		theStream << std::endl;
	}

	return theStream;
}



bool mmpbsa_io::get_next_snap(mmpbsa_io::trajectory_t& traj, std::valarray<mmpbsa_t>& snapshot)
{

#ifdef USE_GROMACS
	if(traj.gromacs_filename != 0)
	{
		mmpbsa_io::load_gmx_trr(*traj.gromacs_filename,snapshot,traj.curr_snap - 1);
		traj.curr_snap++;
		return snapshot.size() != 0;
	}
#endif
	using std::iostream;using std::fstream;
	iostream* sander_file = traj.sander_crd_stream;
	bool returnMe;
	if(sander_file == 0)
	{
		if(traj.sander_filename == 0)
			throw mmpbsa::MMPBSAException("mmpbsa_io::seek: Filename for sander trajectory is a null pointer.",mmpbsa::NULL_POINTER);
		fstream* sander_fstream = new fstream;
		sander_fstream->open(traj.sander_filename->c_str(),std::ios::in);;
		sander_file = sander_fstream;
	}
	sander_file->seekg(traj.curr_pos,sander_file->beg);
	if(traj.natoms == 0 || sander_file == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::get_next_snap: Sander parameters and/or sander coordinate stream is missing.",mmpbsa::DATA_FORMAT_ERROR);
	if(!sander_file->good())
	{
		std::ostringstream error;
		error << "mmpbsa_io::get_next_snap: Cannot obtain " << traj.curr_snap << "th snapshot from trajectory file.";
		throw mmpbsa::MMPBSAException(error,mmpbsa::FILE_IO_ERROR);
	}
	returnMe = get_next_snap(*sander_file,snapshot,traj.natoms,(traj.ifbox > 0));

	traj.curr_pos = sander_file->tellg();
	if(returnMe)
		traj.curr_snap++;
	if(sander_file != traj.sander_crd_stream)
		delete sander_file;

	return returnMe;
}


void mmpbsa_io::seek(mmpbsa_io::trajectory_t& traj,const size_t& snap_pos)
{
#ifdef USE_GROMACS
	if(traj.gromacs_filename != 0)
	{
		traj.curr_snap = snap_pos;
		return;
	}
#endif
	size_t i = traj.curr_snap;
	using std::fstream;using std::iostream;
	iostream* sander_file = traj.sander_crd_stream;
	if(sander_file == 0)
	{
		if(traj.sander_filename == 0)
			throw mmpbsa::MMPBSAException("mmpbsa_io::seek: Filename for sander trajectory is a null pointer.",mmpbsa::NULL_POINTER);
		fstream* sander_fstream = new fstream;
		sander_fstream->open(traj.sander_filename->c_str(),std::ios::in);
		sander_file = sander_fstream;
	}
	sander_file->seekg(traj.curr_pos,sander_file->beg);
	if(traj.natoms == 0)
		throw mmpbsa::MMPBSAException("mmpbsa_io::seek: Trajectory cannot be read without paramters. However, sander paramter object is a null pointer.",mmpbsa::NULL_POINTER);

	for(;i<snap_pos;i++)
	{
		sander_file->seekg(traj.natoms * 24,sander_file->cur);//3 coordinates per atom and 8 characters per coordinate.
		sander_file->seekg(ceil(traj.natoms * 24/80),sander_file->cur);//account for space and end of line characters for each line.
		if(sander_file->fail())
		{
			std::ostringstream error;
			error << "mmpbsa_io::seek: Problem skipping to snap number " << snap_pos;
			throw mmpbsa::MMPBSAException(error,mmpbsa::FILE_IO_ERROR);
		}
	}

	if(traj.ifbox > 0)
		sander_file->seekg(26,sander_file->cur);

	traj.curr_snap = snap_pos;
	traj.curr_pos = sander_file->tellg();
	if(sander_file != traj.sander_crd_stream)
		delete sander_file;

}

void mmpbsa_io::default_trajectory(mmpbsa_io::trajectory_t& traj)
{
	traj.sander_crd_stream = 0;
	traj.sander_filename = 0;
	traj.curr_pos = 0;
	traj.natoms = 0;
	traj.ifbox = 0;

	traj.gromacs_filename = 0;
	traj.curr_snap = 0;
}

void mmpbsa_io::destroy_trajectory(mmpbsa_io::trajectory_t& traj)
{
	delete traj.sander_crd_stream;
	delete traj.sander_filename;
	delete traj.gromacs_filename;
}

mmpbsa_io::trajectory_t mmpbsa_io::open_trajectory(const std::string& filename,const bool& should_remain_in_memory)
{
	using std::fstream;
	trajectory_t returnMe;
	bool is_sander = true;
	mmpbsa_io::default_trajectory(returnMe);

#ifdef USE_GROMACS
	if(filename.find(".trr") != std::string::npos)
	{
		returnMe.gromacs_filename = new std::string(filename);
		return returnMe;
	}
#endif
	//Test whether the file can be opened. If so, setup trajectory. Otherwise, throw exception.
	fstream* trajDiskFile = new fstream(filename.c_str(),std::ios::in);
	if(!trajDiskFile->good())
		throw mmpbsa::MMPBSAException("mmpbsa_io::open_trajectory: Unable to read from trajectory file",mmpbsa::BROKEN_TRAJECTORY_FILE);
	returnMe.sander_filename = new std::string(filename);
	returnMe.curr_snap = 1;
	returnMe.curr_pos = trajDiskFile->tellg();

	if(should_remain_in_memory)
	{
		returnMe.sander_crd_stream = new std::stringstream;
		smart_read(*returnMe.sander_crd_stream,*trajDiskFile,returnMe.sander_filename);
	}
	trajDiskFile->close();
	delete trajDiskFile;
	trajDiskFile = 0;

	return returnMe;

}

bool mmpbsa_io::eof(trajectory_t& traj)
{
	using std::ifstream;
	if(traj.sander_crd_stream == 0)
	{
#ifdef USE_GROMACS
		if(traj.gromacs_filename != 0)
			return mmpbsa_io::gmx_trr_eof(*traj.gromacs_filename,traj.curr_snap);
#endif
		ifstream::streampos eof;
		if(traj.sander_filename == 0)
			throw mmpbsa::MMPBSAException("mmpbsa_io::eof: No trajectory file provided.",mmpbsa::NULL_POINTER);
		ifstream sander_file(traj.sander_filename->c_str());
		if(!sander_file.good())
			return true;
		sander_file.seekg(0,sander_file.end);
		eof = sander_file.tellg();
		sander_file.seekg(traj.curr_pos,sander_file.beg);
		return (traj.curr_pos >= eof);
	}
	return traj.sander_crd_stream->eof();
}


std::string mmpbsa_io::get_traj_title(mmpbsa_io::trajectory_t& traj)
{
	if(traj.sander_crd_stream != 0)
	{
		std::string returnMe;
		traj.curr_pos = 0;
		traj.sander_crd_stream->seekg(traj.curr_pos,std::ios::beg);
		returnMe = mmpbsa_io::get_traj_title(*traj.sander_crd_stream);
		traj.curr_pos = traj.sander_crd_stream->tellg();
		return returnMe;
	}
	else if(traj.sander_filename != 0)
	{
		std::string returnMe;
		std::fstream sander_file(traj.sander_filename->c_str(),std::ios::in);
		returnMe = mmpbsa_io::get_traj_title(sander_file);
		traj.curr_pos = sander_file.tellg();
		return returnMe;
	}
	else if(traj.gromacs_filename != 0)
		return *traj.gromacs_filename;

	throw mmpbsa::MMPBSAException("mmpbsa_io::get_traj_title: no trajectory provided.",mmpbsa::DATA_FORMAT_ERROR);

}

void init(mmpbsa_io::trajectory_t* traj)
{
	mmpbsa_io::default_trajectory(*traj);
}

void destroy(mmpbsa_io::trajectory_t* traj)
{
	mmpbsa_io::destroy_trajectory(*traj);
}

//explicit instantiation
template bool mmpbsa_io::loadValarray<size_t>(std::iostream&, std::valarray<size_t>&,const size_t&, const size_t&, const size_t&);
template bool mmpbsa_io::loadValarray<int>(std::iostream&, std::valarray<int>&,const size_t&, const size_t&, const size_t&);




