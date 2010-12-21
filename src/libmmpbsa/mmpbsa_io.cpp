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

std::string mmpbsa_io::get_traj_title(std::fstream& trajFile)
{
    trajFile.seekg(0,std::ios::beg);
    return getNextLine(trajFile);
}

std::string mmpbsa_io::getNextLine(std::fstream& file) throw (mmpbsa::MMPBSAException)
{
    if(!file.good())
        throw mmpbsa::MMPBSAException("Could not read from file");

    std::string returnMe;
    getline(file,returnMe);
    return returnMe;
}

bool mmpbsa_io::get_next_snap(std::fstream& trajFile, std::valarray<mmpbsa_t>& snapshot,
    const size_t& natoms,bool isPeriodic)
{
    bool returnMe = loadValarray(trajFile,snapshot,natoms*3,8,10);
    if(isPeriodic)
        getNextLine(trajFile);//ignoring periodic box information
    return returnMe;
}

size_t mmpbsa_io::count_snapshots(std::fstream& trajFile,const size_t& natoms, bool isPeriodic)
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

void mmpbsa_io::skip_next_snap(std::fstream& trajFile, const size_t& natoms, bool isPeriodic)
{
    //This might seem excessive just to skip, but I want to verify a snapshot is
    //actually there.
    size_t lineIndex = 0;
    size_t width = 8;//character width of data
    for(size_t dataIndex = 0;dataIndex<natoms*3;)
    {
        if(trajFile.eof())
            throw mmpbsa::SanderIOException("Trajectory file ended in the middle of the "
                    "trajectory.",mmpbsa::UNEXPECTED_EOF);

        std::string currentLine = getNextLine(trajFile);//do not trim string. Spaces are part of formatted size.
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
            dataIndex++;
            currentLine.erase(0,width);
        }

        lineIndex++;
    }
    if(isPeriodic)
        getNextLine(trajFile);
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


template <class T> bool mmpbsa_io::loadValarray(std::fstream& dataFile,
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

template <> bool mmpbsa_io::loadValarray<std::string>(std::fstream& dataFile,
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

void mmpbsa_io::read_siz_file(std::fstream& theFile,
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

int mmpbsa_io::fileopen(const char* filename, const std::ios::openmode& mode, std::fstream& file)
{
    if(file.is_open())
        file.close();
    file.open(filename,mode);
    if(file.is_open())
        return 0;
    else
        return int(mmpbsa::FILE_IO_ERROR);
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

//explicit instantiation
template bool mmpbsa_io::loadValarray<size_t>(std::fstream&, std::valarray<size_t>&,const size_t&, const size_t&, const size_t&);
template bool mmpbsa_io::loadValarray<int>(std::fstream&, std::valarray<int>&,const size_t&, const size_t&, const size_t&);




