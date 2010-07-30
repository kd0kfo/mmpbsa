#include "mmpbsa_io.h"

std::string mmpbsa_io::read_crds(std::fstream& crdFile, std::valarray<mmpbsa_t>& crds)
{
    using std::string;
    using namespace mmpbsa_utils;
    
    if(!crdFile.good())
        throw mmpbsa::SanderIOException("Cannot open coordinate file",mmpbsa::FILE_READ_ERROR);

    string title = getNextLine(crdFile);
    string strNatoms = getNextLine(crdFile);
    strNatoms =trimString(strNatoms);
    size_t natoms = 0;
    sscanf(strNatoms.c_str(),"%d",&natoms);

    if(!loadValarray(crdFile,crds,natoms*3,12,8))
        throw mmpbsa::SanderIOException("Coordinate file is too short.",mmpbsa::FILE_READ_ERROR);

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
    char* format = "%12.7f";//format of the coordinate data

    if(!outFile.good())
    {
        char error[256];
        sprintf(error,"Could not open: %s",fileName);
        throw mmpbsa::SanderIOException(error,mmpbsa::FILE_READ_ERROR);
    }

    char strOutput[12];//used for outputting Fortran formatted strings with sprintf.

    outFile << title << std::endl;

    sprintf(strOutput,"%5d",natoms);//number of atoms
    outFile << strOutput << std::endl;

    size_t m;
    mmpbsa_t dblOutput;
    //save data in rows of 6
    for(m = 0;m<crds.size() - 6;m+=6)
    {
        valarray<mmpbsa_t> row = crds[slice(m,6,1)];//m-th row

        for(size_t i = 0;i<6;i++)
        {
            dblOutput = row[i];
            sprintf(strOutput,format,dblOutput);
            outFile << strOutput;
        }
        outFile << std::endl;
    }

    //save the possibly incomplete last row.
    if(m<crds.size())
    {
        for(m;m<crds.size();m++)
        {
            dblOutput = crds[m];
            sprintf(strOutput,format,dblOutput);
            outFile << strOutput;
        }
        outFile << std::endl;
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
            char* error;
            sprintf(error,"Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #%d is only"
                    "%d characters long.",lineIndex+1,currentLine.size());
            std::cerr << error << std::endl;
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
    mmpbsa_t fltCurrentData = 0;
    size_t dataIndex = 0;

    for(dataIndex;dataIndex<arrayLength;)
    {
        if(dataFile.eof())
            throw mmpbsa::SanderIOException("Data file ended in the middle of the "
                    "data.",mmpbsa::UNEXPECTED_EOF);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
            char* error;
            sprintf(error,"Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #%d is only"
                    "%d characters long.",lineIndex+1,currentLine.size());
            std::cerr << error << std::endl;
        }

        //tokenize line into data. put data into valarray.
        while(currentLine.size() >= width)
        {
            string currentData = currentLine.substr(0,width);
            sscanf(currentData.c_str(),"%f",&fltCurrentData);
            dataArray[dataIndex++] = T(fltCurrentData);
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

    for(dataIndex;dataIndex<arrayLength;)
    {
        if(dataFile.eof())
            throw mmpbsa::SanderIOException("Data file ended in the middle of the "
                    "data.",mmpbsa::BROKEN_TRAJECTORY_FILE);

        string currentLine = getNextLine(dataFile);//do not trim string. Spaces are part of formatted size.
        if(currentLine.size() % width )
        {
            char* error;
            sprintf(error,"Data file contains a short line. "
                    "Lines must be at least 36 characters, but line #%d is only"
                    "%d characters long.",lineIndex+1,currentLine.size());
            std::cerr << error << std::endl;
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

void mmpbsa_io::read_siz_file(std::fstream& theFile,
        std::map<std::string,mmpbsa_t>& radii, std::map<std::string,std::string>& residues)
{
    using mmpbsa_utils::trimString;
    using mmpbsa_utils::toUpperCase;
    
    if(!theFile.good())
        throw mmpbsa::MMPBSAException("Could not open SIZ file.",mmpbsa::FILE_READ_ERROR);

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
            char error[128];
            sprintf(error,"Improperly formatted SIZ file: Short line at %d ",lineNumber);
            throw mmpbsa::MMPBSAException(error,mmpbsa::FILE_READ_ERROR);
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
        data = trimString(currLine.substr(9));
        sscanf(data.c_str(),"%f",&fData);

        radii[atomName] = fData;
        residues[atomName] = residue;
    }

}

int mmpbsa_io::fileopen(const char* filename, const std::ios::openmode& mode, std::fstream& file)
{

#ifdef __USE_BOINC__
    std::string resolved_name;
    int retval = boinc_resolve_filename_s(filename, resolved_name);
    if(retval)
        return retval;

    if(file.is_open())
        file.close();
    file.open(resolved_name.c_str(),mode);
    return retval;
#else
    if(file.is_open())
        file.close();
    file.open(filename,mode);
    if(file.is_open())
        return 0;
    else
        return 1;
#endif

}

//explicit instantiation
template bool mmpbsa_io::loadValarray<size_t>(std::fstream&, std::valarray<size_t>&,const size_t&, const size_t&, const size_t&);
template bool mmpbsa_io::loadValarray<int>(std::fstream&, std::valarray<int>&,const size_t&, const size_t&, const size_t&);




