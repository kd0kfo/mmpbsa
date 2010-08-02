#include "EnergyInfo.h"

mmpbsa::EnergyInfo::EnergyInfo() {
    energydata;
    energydata.resize(total_parameters,0);
}

mmpbsa::EnergyInfo::EnergyInfo(const mmpbsa::EnergyInfo& orig) {
    energydata = orig.energydata;
}

mmpbsa::EnergyInfo::~EnergyInfo() {
}

mmpbsa::EnergyInfo& mmpbsa::EnergyInfo::operator=(const mmpbsa::EnergyInfo& rhs)
{
    if (this == &rhs)
        return *this;

    energydata = rhs.energydata;
    return *this;

}

mmpbsa::EnergyInfo mmpbsa::EnergyInfo::operator+(const mmpbsa::EnergyInfo& rhs)const
{
    EnergyInfo tmp;
    tmp.energydata = energydata+rhs.energydata;
    return tmp;
}

void mmpbsa::EnergyInfo::operator+=(const EnergyInfo& rhs)
{
    energydata += rhs.energydata;
}
mmpbsa::EnergyInfo mmpbsa:: EnergyInfo::operator-(const mmpbsa::EnergyInfo& rhs)const
{
    EnergyInfo tmp;
    tmp.energydata = energydata-rhs.energydata;
    return tmp;
}

void mmpbsa::EnergyInfo::operator-=(const mmpbsa::EnergyInfo& rhs)
{
    energydata -= rhs.energydata;
}

mmpbsa::EnergyInfo mmpbsa::EnergyInfo::operator/(const mmpbsa::EnergyInfo& rhs)const
{
    mmpbsa::EnergyInfo returnMe;
    returnMe.energydata = energydata/rhs.energydata;
    return returnMe;
}

void mmpbsa::EnergyInfo::operator/=(const mmpbsa::EnergyInfo& rhs)
{
    energydata /= rhs.energydata;
}

void mmpbsa::EnergyInfo::get_next_energyinfo(std::fstream& mdoutFile)
{
    using std::string;
    using mmpbsa_io::getNextLine;
    using mmpbsa_utils::trimString;
    
    if(!mdoutFile.good())
        throw SanderIOException("Cannot open mdout file.",FILE_READ_ERROR);

    string currentLine = getNextLine(mdoutFile);
    while(currentLine.size() < 6 || currentLine.substr(1,5) != "NSTEP")//" NSTEP" begins a block of energy info
        currentLine = getNextLine(mdoutFile);
    
    string token;
    while(currentLine.substr(1,5) != "-----")
    {
        currentLine = trimString(currentLine);
        mmpbsa_utils::StringTokenizer tokens(currentLine);
        while(tokens.hasMoreTokens())
        {
            string identifier = tokens.nextToken();
            string tempString = tokens.nextToken();
            while(tempString != "=")
            {
                identifier.append(" ").append(tempString);
                tempString = tokens.nextToken();
            }
            string value = tokens.nextToken();
            if(!loadEnergyValue(identifier,value))
                std::cerr << "Warning: Unknown energy type: " << identifier << std::endl;
        }
        currentLine = getNextLine(mdoutFile);
    }


}

void mmpbsa::EnergyInfo::get_first_energyinfo(const char* fileName)
{
    std::fstream mdout;
    mmpbsa_io::fileopen(fileName,std::ios::in,mdout);
    get_next_energyinfo(mdout);
    mdout.close();
}

bool mmpbsa::EnergyInfo::loadEnergyValue(const std::string& identifier,const std::string& value)
{
    float dblValue = 0;
    sscanf(value.c_str(),MMPBSA_FORMAT,&dblValue);
    return loadEnergyValue(identifier,dblValue);
}

bool mmpbsa::EnergyInfo::loadEnergyValue(const std::string& identifier,const mmpbsa_t& value)
{
    //ensure the energydata valarray is the correct size. If not, whatever is
    //there is lost. Not a problem, because in that case something is wrong anyways.
    if(energydata.size() != this->total_parameters)
        energydata.resize(total_parameters,0);

    //Long if..else chain ahead. Sorry, have a lot to check. No string switch in c++. 
    //Well, one could make one with map, but still.
    if(identifier == "NSTEP")
        energydata[nstep] = value;
    else if(identifier == "TIME(PS)")
        energydata[time] = value;
    else if(identifier == "TEMP(K)")
        energydata[temp] = value;
    else if(identifier == "PRESS")
        energydata[press] = value;
    else if(identifier == "Etot")
        energydata[etot] = value;
    else if(identifier == "EKtot")
        energydata[ektot] = value;
    else if(identifier == "EPtot")
        energydata[eptot] = value;
    else if(identifier == "BOND")
        energydata[bond] = value;
    else if(identifier == "ANGLE")
        energydata[angle] = value;
    else if(identifier == "DIHED")
        energydata[dihed] = value;
    else if(identifier == "1-4 NB")
        energydata[nb14] = value;
    else if(identifier == "1-4 EEL")
       energydata[eel14] = value;
    else if(identifier == "VDWAALS")
       energydata[vdwaals] = value;
    else if(identifier == "EELEC")
        energydata[eelec] = value;
    else if(identifier == "EHBOND")
        energydata[ehbond] = value;
    else if(identifier == "EGB")
        energydata[egb] = value;
    else if(identifier == "RESTRAINT")
        energydata[restraint] = value;
    else if(identifier == "ESURF")
        energydata[esurf] = value;
    else if(identifier == "EAMBER (non-restraint)")
        energydata[nonconst_pot] = value;
    else if(identifier == "DV/DL")
        energydata[dvdl] = value;
    else if(identifier == "|E(PBS)")
        energydata[rms_pbs] = value;
    else if(identifier == "EKCMT")
        energydata[ekcmt] = value;
    else if(identifier == "VIRIAL")
        energydata[virial] = value;
    else if(identifier == "VOLUME")
        energydata[volume] = value;
    else if(identifier == "T_non-LES")
        energydata[eksolt] = value;
    else if(identifier == "T_LES")
        energydata[eksolv] = value;
    else if(identifier == "EPOLZ")
        energydata[epol] = value;
    else if(identifier == "E3BODY")
        energydata[e3bod] = value;
    else if(identifier == "Dipole convergence: rms")
        energydata[diprms] = value;
    else if(identifier == "iters")
        energydata[dipitr] = value;
    else if(identifier == "temperature")
        energydata[dipole_temp] = value;
    else if(identifier == "Density")
        energydata[density] = value;
    else if(identifier == "Ewald error estimate")
        energydata[ewalderr] = value;
    else if(identifier == "Current RMSD from reference")
        energydata[rmsdvalue] = value;
    else if(identifier == "Current target RMSD" )
        energydata[tgtrmsd] = value;
    else
        return false;

    return true;
}

void mmpbsa::EnergyInfo::clear()
{
    energydata *= 0.0;
}

mmpbsa::AveRmsEnerInfo::AveRmsEnerInfo()
{
    avg;
    rms;
    relrms;
}

mmpbsa::AveRmsEnerInfo::AveRmsEnerInfo(const mmpbsa::AveRmsEnerInfo& orig)
{
    avg = orig.avg;
    rms = orig.rms;
    relrms = orig.relrms;
}

mmpbsa::AveRmsEnerInfo::AveRmsEnerInfo(const mmpbsa::EnergyInfo& avgs, const mmpbsa::EnergyInfo& rmses)
{
    using std::valarray;

    //Directly copy avg and rms
    avg = avgs;
    rms = rmses;
    relrms;

    //Calculate relrms, according to relrms = rms/abs(avg)
    const valarray<mmpbsa_t>& avgdata = avg.getEnergyData();
    const valarray<mmpbsa_t>& rmsdata = rms.getEnergyData();
    valarray<mmpbsa_t> relrmsData(0.0,avgdata.size());
    relrmsData[0] = rmsdata[0];
    for(int i = 1;i<avgdata.size();i++)
        if(avgdata[i] != 0)
            relrmsData[i] = rmsdata[i]/std::abs(avgdata[i]);
}

mmpbsa::AveRmsEnerInfo& mmpbsa::AveRmsEnerInfo::operator=(const mmpbsa::AveRmsEnerInfo& rhs)
{
    if(this == &rhs)
        return *this;

    avg = rhs.avg;
    rms = rhs.rms;
    relrms = rhs.relrms;
    return *this;
}

void mmpbsa::AveRmsEnerInfo::get_first_energyinfo(const char* fileName)
{
    1;
}

    /**
 * Loads next energy info from the provided mdout file.
 *
 * @param mdoutFile
 * @return
 */
void mmpbsa::AveRmsEnerInfo::get_next_energyinfo(std::fstream& mdoutFile)
{
    using std::string;
    using mmpbsa_io::getNextLine;
    
    string currentLine = getNextLine(mdoutFile);
    while(currentLine.find("A V E R A G E",0) == currentLine.npos)// "A V E R A G E" begins a block of energy info
        currentLine = getNextLine(mdoutFile);

    EnergyInfo avgs,rmses;
    avgs.get_next_energyinfo(mdoutFile);
    rmses.get_next_energyinfo(mdoutFile);
    *this = AveRmsEnerInfo(avgs,rmses);
}
void mmpbsa::AveRmsEnerInfo::clear()
{
    avg.clear();
    rms.clear();
    relrms.clear();
}

void mmpbsa::mdout2enerinfos(std::fstream& mdout, std::valarray<mmpbsa::EnergyInfo>& energyinfos,
        std::valarray<mmpbsa::AveRmsEnerInfo>& avginfos)
{
    /**
    """Given an open file object for an MDOUT file, return a pair of
    ListDicts corresponding to all the ``energy information'' entries
    in the MDOUT file of an MD run.  The first ListDict contains the
    info from the 'normal' entries written every NTPR steps.  The
    second listdict contains any Ave and RMS info (as in the
    AveRmsEnerInfo class) found in the file.
    """
     */

    using std::string;
    using mmpbsa_io::getNextLine;

    //For valarray, the size needs to be known. This section scans the file to
    //calculate the number of EnergyInfo's and AveRmsEnerInfo's that will be needed.
    //Reading the file twice bothers me, but I'm not sure which is worse,
    //reading the file twice, expanding the valarray everytime, or using vector
    //and then shlepping everything over to the valarray (which I want for arithmetic).
    mdout.seekg(0,std::ios::beg);
    string currentLine = "";
    int numberOfEnergyInfos = 0;
    int numberOfAvgEnergyInfos = 0;
    while(!mdout.eof())
    {
        currentLine = getNextLine(mdout);
        if(currentLine.size() < 6 && currentLine.substr(1,5) != "NSTEP")
            numberOfEnergyInfos++;
        else if(currentLine.find("A V E R A G E",0) != currentLine.npos)
            numberOfAvgEnergyInfos++;
    }
    mdout.seekg(0,std::ios::beg);

    //set proper valarray size
    if(energyinfos.size() != numberOfEnergyInfos)
        energyinfos.resize(numberOfEnergyInfos);
    if(avginfos.size() != numberOfAvgEnergyInfos)
        avginfos.resize(numberOfAvgEnergyInfos);

    //read file and fill arrays
    int avginfoIndex = 0;
    int energyinfoIndex = 0;

    mmpbsa::EnergyInfo newEnergyInfo;
    mmpbsa::AveRmsEnerInfo newAvgInfo;
    while(!mdout.eof())
    {
        //peek
        currentLine = getNextLine(mdout);
        mdout.seekg(-1*currentLine.size(),std::ios::cur);

        //decide what to do
        if(currentLine.size() < 6 && currentLine.substr(1,5) != "NSTEP")
        {
            newEnergyInfo.clear();
            newEnergyInfo.get_next_energyinfo(mdout);
            energyinfos[energyinfoIndex++] = newEnergyInfo;
        }
        else if(currentLine.find("A V E R A G E",0) != currentLine.npos)
        {
            newAvgInfo.clear();
            newAvgInfo.get_next_energyinfo(mdout);
            avginfos[avginfoIndex++] = newAvgInfo;
        }
        else
            getNextLine(mdout);//advance to next line to skip other information
    }
            
}

float mmpbsa::get_minimized_energy(std::fstream& mdout) throw (SanderIOException)
{
    using std::string;
    using mmpbsa_io::getNextLine;
    
    if(!mdout.is_open())
        throw mmpbsa::SanderIOException("Cannot read mdout file to obtain final energy",FILE_READ_ERROR);

    mdout.seekg(0,std::ios::beg);

    //Find final result block. Throw exception if it doesn't exist.
    string currentLine = getNextLine(mdout);
    while(currentLine.find("FINAL RESULTS") == currentLine.npos)
    {
        if(mdout.eof())
            throw SanderIOException("FINAL RESULTS block was missing from the file",FILE_READ_ERROR);
        currentLine = getNextLine(mdout);
    }

    while(currentLine.find("ENERGY") == currentLine.npos)
    {
        if (mdout.eof())
            throw SanderIOException("Total Energy is missing from the FINAL RESULTS block.", FILE_READ_ERROR);
        currentLine = getNextLine(mdout);
    }
    currentLine = getNextLine(mdout);
    mmpbsa_utils::StringTokenizer tokens(currentLine);
    tokens.nextToken();//NSTEP
    string strEnergy = tokens.nextToken();
    float fEnergy = 0;
    sscanf(strEnergy.c_str(),MMPBSA_FORMAT,&fEnergy);

    return fEnergy;

}

void mmpbsa::EnergyInfo::setEnergyData(const std::valarray<mmpbsa_t>& newData)
{
    if(newData.size() != this->total_parameters)
    {
        char error[256];
        sprintf(error,"EnergyInfo needs %d data values, but %d were provided.",
                this->total_parameters,newData.size());
        throw MMPBSAException(error,INVALID_ARRAY_SIZE);
    }
    if(energydata.size() != this->total_parameters)
        energydata.resize(this->total_parameters);
    energydata = newData;
}

