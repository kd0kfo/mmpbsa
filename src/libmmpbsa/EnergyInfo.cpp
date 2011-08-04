#include <iomanip>

#include <fstream>
#include <iostream>

#include "mmpbsa_io.h"
#include "mmpbsa_utils.h"
#include "mmpbsa_exceptions.h"
#include "StringTokenizer.h"
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

void mmpbsa::EnergyInfo::clear(){*this *= 0.0;}

mmpbsa::EnergyInfo sqrt(const mmpbsa::EnergyInfo& rhs)
{
  mmpbsa::EnergyInfo returnMe = rhs;
  for(size_t i = 0;i<returnMe.size();i++)
    returnMe[i] = sqrt(returnMe[i]);
  return returnMe;
}

mmpbsa::EnergyInfo abs(const mmpbsa::EnergyInfo& rhs)
{
  mmpbsa::EnergyInfo returnMe = rhs;
  for(size_t i = 0;i<returnMe.size();i++)
    returnMe[i] = std::abs(returnMe[i]);
  return returnMe;
}

std::fstream& operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data)
{
  int retval = data.get_next_energyinfo(theStream);
  if(retval == mmpbsa::UNEXPECTED_EOF)
    throw mmpbsa::SanderIOException("operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data) reached an unexpected end of file.",(mmpbsa::MMPBSAErrorTypes)retval);
  else if(retval > 0 && retval < mmpbsa::NUMBER_OF_ERROR_TYPES)
    throw mmpbsa::SanderIOException("operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data) encountered an error.",(mmpbsa::MMPBSAErrorTypes)retval);
  else if(retval != 0)
    throw mmpbsa::MMPBSAException("operator>>(std::fstream& theStream, mmpbsa::EnergyInfo& data) encountered an unknown error.");
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
  theStream << "1-4 VDW =  " << data[EnergyInfo::nb14];
  theStream << "  1-4 EEL =   " << data[EnergyInfo::eel14];
  theStream << "  VDWAALS    =  " << data[EnergyInfo::vdwaals];
  theStream << std::endl;
  theStream << "EEL  =  "  << data[EnergyInfo::eelec];
  theStream << "  HBOND  =   "  << data[EnergyInfo::eelec];
  theStream << "  RESTRAINT  =  "  << data[EnergyInfo::restraint];
  theStream << std::endl;
  theStream << "EKCMT  =  " << data[EnergyInfo::ekcmt];
  theStream << "  VIRIAL  =   " << data[EnergyInfo::virial];
  theStream << "  VOLUME     =    " << data[EnergyInfo::volume];
  theStream << std::endl;
  theStream << "Density    =    " << data[EnergyInfo::density];
  theStream << std::endl;
  theStream << "EAMBER (non-restraint)  =  " << data[EnergyInfo::nonconst_pot];
  theStream << std::endl;
  theStream << "Ewald error estimate:   " << data[EnergyInfo::ewalderr];

  return theStream;
}

mmpbsa_t* mmpbsa::EnergyInfo::get_minimization_header(const std::string& header_line) throw (mmpbsa::SanderIOException)
{
  std::istringstream buff;
  mmpbsa_t* returnMe,temp_var;
  size_t delim_loc;
  std::string error_message,_the_line = mmpbsa_utils::trimString(header_line);
  if(_the_line.size() == 0)
    throw mmpbsa::SanderIOException("mmpbsa::EnergyInfo::get_minimization_energy: the line which should contain"
				    " energy information for a minimization contains no data values: " + header_line,
				    mmpbsa::DATA_FORMAT_ERROR);


  error_message = "mmpbsa::EnergyInfo::get_minimization_energy: the line which should contain"
    " energy information for a minimization contains unexpected values: " + _the_line;
  returnMe = new mmpbsa_t[2];
  for(int i = 0;i<2;i++)
    {
      delim_loc = _the_line.find(" ",0);
      if(delim_loc == std::string::npos)
	{
	  delete [] returnMe;
	  throw mmpbsa::SanderIOException(error_message + "\nDelimiter not Found",mmpbsa::DATA_FORMAT_ERROR);
	}
      buff.clear();buff.str(_the_line.substr(0,delim_loc));
      _the_line.erase(0,delim_loc);
      _the_line = mmpbsa_utils::trimString(_the_line);
      buff >> std::scientific >> temp_var;
      if(buff.fail())
	throw mmpbsa::SanderIOException(error_message + "\nData format error: " + buff.str() + "",mmpbsa::DATA_FORMAT_ERROR);;
      returnMe[i] = temp_var;
    }
  return returnMe;
}

int mmpbsa::EnergyInfo::get_next_energyinfo(std::fstream& mdoutFile)
{
  using std::string;
  using mmpbsa_io::getNextLine;
  using mmpbsa_utils::trimString;
   
  size_t line_count = 0;
  if(!mdoutFile.good())
    throw SanderIOException("mmpbsa::EnergyInfo::get_next_energyinfo: Cannot open mdout file.",FILE_IO_ERROR);

  string currentLine = getNextLine(mdoutFile);line_count++;
  while(!mdoutFile.eof() && currentLine.find("NSTEP") == std::string::npos)//" NSTEP" begins a block of energy info
    {
      if(currentLine.find(MDOUT_AVERAGE_HEADER) != string::npos)
	return mmpbsa::UNEXPECTED_EOF;
      currentLine = getNextLine(mdoutFile);line_count++;
    }
    
  if(mdoutFile.eof())
    return mmpbsa::UNEXPECTED_EOF;

  //The format of MDOUT varies depending on whether or not the run was minimization :-/
  //Therefore, we must parse the NSTEP line differently if it is a minimization.
  if(currentLine.find("=") == string::npos)//in this case, this is minimization
    {
      mmpbsa_t* minimization_header;
      currentLine = getNextLine(mdoutFile);line_count++;
      try{
	minimization_header = get_minimization_header(currentLine);
      }
      catch(mmpbsa::SanderIOException sioe)
    	{
	  std::ostringstream error;
	  error << "mmpbsa::EnergyInfo::get_next_energyinfo: Error with minimization data."
		<< std::endl << sioe;
	  throw mmpbsa::SanderIOException(error,sioe.getErrType());
    	}
      (*this)[EnergyInfo::nstep] = minimization_header[0];
      (*this)[EnergyInfo::etot] = minimization_header[1];
      delete [] minimization_header;
      getNextLine(mdoutFile);//minimization skips a line before the Energy = ... block.

    }

  currentLine = getNextLine(mdoutFile);line_count++;
  //NSTEP data format.
  //Line with titles of information
  //Line with information about the STEP
  //Blank Line
  //Energy values of the form "NAME = DATA"
  while(!mdoutFile.eof())
    {
      currentLine = trimString(currentLine);
      if(currentLine.size() == 0 || currentLine.find("------------------------------------") != std::string::npos)
	break;//end of data section.
      if(currentLine.find(CR_CHAR) != std::string::npos)
        {
	  currentLine.erase(currentLine.find(CR_CHAR),1);
	  if(currentLine.size() == 0)
	    break;
        }
      mmpbsa_utils::StringTokenizer tokens(currentLine);
      try{
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
      }
      catch(mmpbsa::TokenizerException te)
	{
	  std::ostringstream error;
	  error << te.what() << std::endl << "On line " << line_count << ": " << currentLine;
	  throw mmpbsa::TokenizerException(error.str(),te.getErrType());
	}

      if(!mdoutFile.eof())
	{
	  currentLine = getNextLine(mdoutFile);line_count++;
	}
      else
	break;
    }

  return 0;

}

int mmpbsa::EnergyInfo::get_first_energyinfo(const char* fileName)
{
  std::fstream mdout(fileName,std::ios::in);
  int retval = get_next_energyinfo(mdout);
  mdout.close();
  return retval;
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
  else if(identifier == "1-4 VDW" || identifier == "1-4 NB")
    energydata[nb14] = value;
  else if(identifier == "1-4 EEL")
    energydata[eel14] = value;
  else if(identifier == "VDWAALS")
    energydata[vdwaals] = value;
  else if(identifier == "EEL" || identifier == "EELEC")
    energydata[eelec] = value;
  else if(identifier == "HBOND" || identifier == "EHBOND")
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
int mmpbsa::AveRmsEnerInfo::get_avg_rms_info(std::fstream& mdoutFile)
{
  using std::string;
  using mmpbsa_io::getNextLine;
    
  string currentLine = getNextLine(mdoutFile);
  EnergyInfo avgs,rmses;
  int avg_retval, rmses_retval;
  while(currentLine.find(MDOUT_AVERAGE_HEADER,0) == currentLine.npos)// "A V E R A G E" begins a block of energy info
    currentLine = getNextLine(mdoutFile);

  avg_retval = avgs.get_next_energyinfo(mdoutFile);
  rmses_retval = rmses.get_next_energyinfo(mdoutFile);
  *this = AveRmsEnerInfo(avgs,rmses);

  return (avg_retval != 0) ? avg_retval : rmses_retval;
}
void mmpbsa::AveRmsEnerInfo::clear()
{
  avg.clear();
  rms.clear();
  relrms.clear();
}

mmpbsa_t mmpbsa::get_minimized_energy(std::fstream& mdout) throw (SanderIOException)
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

const mmpbsa::EnergyInfo& mmpbsa::AveRmsEnerInfo::get_average()const{return avg;}
const mmpbsa::EnergyInfo& mmpbsa::AveRmsEnerInfo::get_rms()const{return rms;}
const mmpbsa::EnergyInfo& mmpbsa::AveRmsEnerInfo::get_relrms()const{return relrms;}


