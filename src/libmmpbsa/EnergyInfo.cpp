#include "EnergyInfo.h"

bool mmpbsa::EnergyInfo::operator>(const mmpbsa::EnergyInfo& rhs)const
{
    if(rhs.size() != this->size())
        return false;
    for(size_t i = 0;i<this->size();i++)
        if((*this)[i] > rhs[i])
            return true;

    return false;
}

mmpbsa::EnergyInfo sqrt(const mmpbsa::EnergyInfo& rhs)
{
    mmpbsa::EnergyInfo returnMe = rhs;
    for(size_t i = 0;i<returnMe.size();i++)
        returnMe[i] = sqrt(returnMe[i]);
    return returnMe;
}

std::fstream& operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data)
{
    data.get_next_energyinfo(theStream);
    return theStream;
}

std::ostream& operator<<(std::ostream& theStream, const mmpbsa::EnergyInfo& data)
{
	using mmpbsa::EnergyInfo;
    if(data.size() != EnergyInfo::total_parameters)
        throw mmpbsa::MMPBSAException("Uninitialized EnergyInfo object in operator<<.",mmpbsa::DATA_FORMAT_ERROR);
    
    theStream << "NSTEP = " <<  data[EnergyInfo::nstep];
    theStream << "  TIME(PS) = " << data[EnergyInfo::time];
    theStream << "  TEMP(K) =  " << data[EnergyInfo::temp];
    theStream << "  PRESS =    " <<data[EnergyInfo::press];
    theStream << std::endl;
    theStream << "Etot   = " << data[EnergyInfo::etot];
    theStream << "  EKtot   =  " << data[EnergyInfo::ektot];
    theStream << "  EPtot   =  " << data[EnergyInfo::eptot];
    theStream << std::endl;
    theStream << "BOND   =  " << data[EnergyInfo::bond];
    theStream << "  ANGLE   =   " << data[EnergyInfo::angle];
    theStream << "  DIHED   =   " << data[EnergyInfo::dihed];
    theStream << std::endl;
    theStream << "1-4 NB =  " << data[EnergyInfo::nb14];
    theStream << "  1-4 EEL =   " << data[EnergyInfo::eel14];
    theStream << "  VDWAALS    =  " << data[EnergyInfo::vdwaals];
    theStream << std::endl;
    theStream << "EELEC  =  "  << data[EnergyInfo::eelec];
    theStream << "  EHBOND  =   "  << data[EnergyInfo::eelec];
    theStream << "  RESTRAINT  =  "  << data[EnergyInfo::restraint];
    theStream << std::endl;
    theStream << "EAMBER (non-restraint)  =  " << data[EnergyInfo::nonconst_pot];
    theStream << std::endl;
    theStream << "Ewald error estimate:   " << data[EnergyInfo::ewalderr];

    return theStream;
}

void mmpbsa::EnergyInfo::get_next_energyinfo(std::fstream& mdoutFile)
{
    using std::string;
    using mmpbsa_io::getNextLine;
    using mmpbsa_utils::trimString;
    
    if(!mdoutFile.good())
        throw SanderIOException("Cannot open mdout file.",FILE_IO_ERROR);

    string currentLine = getNextLine(mdoutFile);
    while(!mdoutFile.eof() && (currentLine.size() < 6 || currentLine.substr(1,5) != "NSTEP"))//" NSTEP" begins a block of energy info
        currentLine = getNextLine(mdoutFile);
    
    string token;
    while(!mdoutFile.eof() && currentLine.substr(1,5) != "-----")
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
                if(*(identifier.end()-1) == ':')
                {
                    identifier.erase(identifier.end()-1);
                    tempString = "=";
                }
                else
                    tempString = tokens.nextToken();
            }
            string value = tokens.nextToken();
            if(!loadEnergyValue(identifier,value))
                std::cerr << "Warning: Unknown energy type: " << identifier << std::endl;
        }
        if(!mdoutFile.eof())
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
    mmpbsa_t dblValue = 0;
    std::istringstream buff(value);
    buff >> MMPBSA_FORMAT >> dblValue;
    return loadEnergyValue(identifier,dblValue);
}

bool mmpbsa::EnergyInfo::loadEnergyValue(const std::string& identifier,const mmpbsa_t& value)
{
    EnergyInfo& energydata = *this;
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

mmpbsa::AveRmsEnerInfo::AveRmsEnerInfo()
{
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
    relrms.resize(avgs.size(),0.0);

    //Calculate relrms, according to relrms = rms/abs(avg)
    relrms[0] = rms[0];
    for(size_t i = 1;i<avg.size();i++)
        if(avg[i] != 0)
            relrms[i] = rms[i]/std::abs(avg[i]);
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
    avg *= 0.0;
    rms *= 0.0;
    relrms *= 0.0;
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
    size_t numberOfEnergyInfos = 0;
    size_t numberOfAvgEnergyInfos = 0;
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
            newEnergyInfo *= 0.0;
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
        throw mmpbsa::SanderIOException("Cannot read mdout file to obtain final energy",FILE_IO_ERROR);

    mdout.seekg(0,std::ios::beg);

    //Find final result block. Throw exception if it doesn't exist.
    string currentLine = getNextLine(mdout);
    while(currentLine.find("FINAL RESULTS") == currentLine.npos)
    {
        if(mdout.eof())
            throw SanderIOException("FINAL RESULTS block was missing from the file",FILE_IO_ERROR);
        currentLine = getNextLine(mdout);
    }

    while(currentLine.find("ENERGY") == currentLine.npos)
    {
        if (mdout.eof())
            throw SanderIOException("Total Energy is missing from the FINAL RESULTS block.", FILE_IO_ERROR);
        currentLine = getNextLine(mdout);
    }
    currentLine = getNextLine(mdout);
    mmpbsa_utils::StringTokenizer tokens(currentLine);
    tokens.nextToken();//NSTEP
    string strEnergy = tokens.nextToken();
    mmpbsa_t energy = 0;
    std::istringstream buff(strEnergy);
    buff >> MMPBSA_FORMAT >> energy;

    return energy;

}



