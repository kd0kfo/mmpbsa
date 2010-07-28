#include "mmpbsa.h"


int main(int argc, char** argv)
{
    try 
    {
#ifdef __BOINC__
        boinc_init();//must be called before any other BOINC routines.
#endif
        currState;
        int returnMe = realDeal(argc, argv);
        currState.flush();
        currState.close();
        return returnMe;
    }    
    catch (mmpbsa::MMPBSAException e)
    {
        std::cerr << e.identifier() << ": " << e.what() << std::endl;
        if(currState.trustPrmtop)
        {
            std::cerr << "The trust_prmtop flag was set. Perhaps that is a problem." << std::endl;
        }
        return e.getErrType();
    }
}

int realDeal(int argc, char** argv)
{
    using std::valarray;
    using std::vector;
    using std::slice;
    using std::map;
    using namespace mmpbsa;

    MeadInterface mi;

    if(argc == 1)
      {
	parseFlag("help",mi);
	return 0;
      }

    currState.receptorStartPos.push_back(0);//would be better as a command line argument, but I don't know yet what this will look like. To add later.
    currState.ligandStartPos.push_back(1);//array containing the starting positions of ligands and receiptors
    
    int retval = parseArgs(argc,argv,mi);
    if(retval)
        return retval - 1;//this way "help, for example, will trigger exiting (retval = 1), but will return the program with an exit value of 0 (i.e. not an error-based exit);

    //load and check the parmtop file.
    mmpbsa::SanderParm * sp = new mmpbsa::SanderParm;
    sp->raw_read_amber_parm(currState.prmtopFile);
    if(!currState.trustPrmtop)
        if(!sp->sanityCheck())
            throw MMPBSAException("Parmtop file is insane.",INVALID_PRMTOP_DATA);
    currState.prmtopFile.close();

    //Create energy function with the parmtop data. This energy function will
    //have everything in it. Receptor and ligand will be stripped out.
    EmpEnerFun entireEFun(sp);

    valarray<bool> complexKeepers(false,sp->natom);//array of atoms to keep.
    valarray<bool> receptorKeepers(false,sp->natom);
    valarray<bool> ligandKeepers(false,sp->natom);
    size_t bottom,top;
    size_t receptorSize = 0;
    size_t ligandSize = 0;

    for(size_t i = 0;i<currState.receptorStartPos.size();i++)
    {
        size_t currPos = currState.receptorStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currReceptor(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorSize += top-bottom;
    }
    for(size_t i = 0;i<currState.ligandStartPos.size();i++)
    {
        size_t currPos = currState.ligandStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currLigand(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currLigand;
        ligandKeepers[slice(bottom,top-bottom,1)] = currLigand;
        ligandSize += top-bottom;
    }
    size_t complexSize = receptorSize+ligandSize;

    EmpEnerFun complexEFun = entireEFun.stripEnerFun(complexKeepers,true);
    EmpEnerFun receptorEFun = entireEFun.stripEnerFun(receptorKeepers,true);
    EmpEnerFun ligandEFun = entireEFun.stripEnerFun(ligandKeepers,true);

    //load radii data, if available
    map<std::string,mmpbsa_t> radii;//later, check to see if radii.size() > 0 before calling full_EMap(...)
    map<std::string,std::string> residues;
    if(currState.radiiFile.is_open())
    {
        mmpbsa_io::read_siz_file(currState.radiiFile,radii, residues);
        currState.radiiFile.close();
    }

    if(!currState.trajFile.good())
        throw MMPBSAException("Unable to read from trajectory file",BROKEN_TRAJECTORY_FILE);

    using namespace mmpbsa_io;
    get_traj_title(currState.trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
    valarray<mmpbsa_t> snapshot(sp->natom*3);
    valarray<mmpbsa_t> complexSnap(complexSize*3);
    valarray<mmpbsa_t> receptorSnap(receptorSize*3);
    valarray<mmpbsa_t> ligandSnap(ligandSize*3);

    bool isPeriodic = sp->ifbox > 0;//Are periodic boundary conditions used?

    //if the program is resuming a previously started calculation, advance to the
    //last snapshot.
    if(currState.currentSnap)//zero = one = start from beginning.
        for(size_t i = 0;i<currState.currentSnap-1;i++)
            try
            {
                mmpbsa_io::skip_next_snap(currState.trajFile,sp->natom,isPeriodic);
            }
            catch(MMPBSAException e)
            {
                if(e.getErrType() == UNEXPECTED_EOF)
                {
                  fprintf(stdout,"End of Snapshots Reached\n");
                  return 0;
                }
            }
    size_t snapcounter = (currState.currentSnap) ? currState.currentSnap - 1 : 0;//snapcounter will be incremented below

    //Walk through the snapshots. This is where MMPBSA is actually done.
    while(!currState.trajFile.eof())
    {
        try{
            //if a list of snaps to be run is provided, check to see if this snapshot
            //should be used. Remember: snapcounter is 1-indexed.
            snapcounter++;
            if(currState.snapList.size())//check if the current snapshot should be skipped
                if(!mmpbsa_utils::contains(currState.snapList,snapcounter))
                {
                    mmpbsa_io::skip_next_snap(currState.trajFile,sp->natom,isPeriodic);
                    printf("Skipping Snapshot #%d\n",snapcounter);
                    continue;
                }
        
            if(get_next_snap(currState.trajFile, snapshot, sp->natom,isPeriodic))
                printf("Running Snapshot #%d\n",snapcounter);
            else
            {
                char error[256];
                sprintf(error,"Error in loading snapshot #%d",++snapcounter);
                throw MMPBSAException(error,BROKEN_TRAJECTORY_FILE);
            }
        }
        catch(MMPBSAException e)
        {
            if(e.getErrType() == UNEXPECTED_EOF)
            {
              std::cerr << "End of Snapshots Reached" << std::endl;
              return 0;
            }
        }

        //process snapshot.
        size_t complexCoordIndex = 0;
        size_t receptorCoordIndex = 0;
        size_t ligandCoordIndex = 0;
        for(size_t i = 0;i<sp->natom;i++)
        {
            std::slice_array<mmpbsa_t> currCoord = snapshot[slice(3*i,3,1)];
            if(complexKeepers[i])
                complexSnap[slice(3*complexCoordIndex++,3,1)] = currCoord;
            if(receptorKeepers[i])
                receptorSnap[slice(3*receptorCoordIndex++,3,1)] = currCoord;
            if(ligandKeepers[i])
                ligandSnap[slice(3*ligandCoordIndex++,3,1)] = currCoord;
        }

        FinDiffMethod fdm = MeadInterface::createFDM(complexSnap,receptorSnap,ligandSnap);
        const map<std::string,mmpbsa_t>* pradii = &(mi.brad);//don't delete!!!
        if(radii.size())//if the radius map is empty, use MeadInterface's lookup table.
            pradii = &radii;

        //output-ing is broken up by section, in case the program needs to be
        //monitored or paused.
        if(currState.currentMolecule == MMPBSAState::END_OF_MOLECULES)
        {
            currState.currentMolecule = MMPBSAState::COMPLEX;
            continue;//Restarted program at the end of a snapshot. So, move on.
        }
        
        if(currState.currentMolecule == MMPBSAState::COMPLEX)
        {
            EMap complexEMap = MeadInterface::full_EMap(complexEFun,complexSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            currState.outputFile << "COMPLEX" << endl << complexEMap << endl;
            currState.currentMolecule = MMPBSAState::RECEPTOR;
            currState.currentSnap = snapcounter;
            currState.checkpoint_out();
        }
        

        if(currState.currentMolecule == MMPBSAState::RECEPTOR)
        {
            EMap receptorEMap = MeadInterface::full_EMap(receptorEFun,receptorSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            currState.outputFile << "RECEPTOR" << endl << receptorEMap << endl;
            currState.currentMolecule = MMPBSAState::LIGAND;
            currState.currentSnap = snapcounter;
            currState.checkpoint_out();
        }

        if(currState.currentMolecule == MMPBSAState::LIGAND)
        {
            EMap ligandEMap = MeadInterface::full_EMap(ligandEFun,ligandSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            currState.outputFile << "LIGAND" << endl << ligandEMap << endl << endl;
            currState.currentMolecule = MMPBSAState::END_OF_MOLECULES;
        }
        currState.checkpoint_out();
        currState.currentMolecule = MMPBSAState::COMPLEX;//Rest current molecule
    }//end of snapshot loop


    delete sp;
    return 0;
}

int parseArgs(int argc, char** argv, mmpbsa::MeadInterface& mi)
{
    using std::string;
    for(int i = 1;i<argc;i++)
    {
        string currArg = argv[i];
        int retval = 0;
        if(currArg.substr(0,2) == "--")
            currArg.erase(currArg.begin(),currArg.begin()+2);

        if(currArg.find("=") != string::npos)
        {
            retval = parseParameter(currArg,mi);
            if(retval)
                return retval;
        }
        else
        {
            retval = parseFlag(currArg,mi);
            if(retval)
                return retval;
        }
    }
    return 0;
}

int parseParameter(std::string arg, mmpbsa::MeadInterface& mi)
{
    if(arg.find("=") == string::npos)
        return parseFlag(arg,mi);//oops that's a flag.
    
    if(arg.find("=") != arg.rfind("="))
    {
        char error[256];
        sprintf(error,"Multiple occurance of \"=\" in parameter: %s\n",arg.c_str());
        throw mmpbsa::MMPBSAException(error,mmpbsa::COMMAND_LINE_ERROR);
    }
    
    using std::string;
    using mmpbsa_io::fileopen;
    string name = arg.substr(0,arg.find("="));
    string value = arg.substr(arg.find("=")+1);

    if(name == "prmtop" || name == "parmtop")
    {
        fileopen(value.c_str(),std::ios::in,currState.prmtopFile);
    }
    else if(name == "coordinates" || name == "traj" || name == "mdcrd" || name == "mdcrds")
    {
        fileopen(value.c_str(),std::ios::in,currState.trajFile);
    }
    else if(name == "out" || name == "output")
    {
        //will append for the sake of restarting in the middle of calculations.
        fileopen(value.c_str(),std::ios::out | std::ios::app,currState.outputFile);
    }
    else if(name == "radii")
    {
        fileopen(value.c_str(),std::ios::in,currState.radiiFile);
    }
    else if(name == "istrength")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.istrength = mmpbsa_t(fValue);
    }
    else if(name == "surf_offset")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.surf_offset = mmpbsa_t(fValue);
    }
    else if(name == "surf_tension")
    {
        float fValue = 0;
        sscanf(value.c_str(),"%f",&fValue);
        mi.surf_tension = mmpbsa_t(fValue);
    }
    else if(name == "rec_list")
    {
        currState.receptorStartPos.clear();
        loadListArg(value,currState.receptorStartPos);
    }
    else if(name == "lig_list")
    {
        currState.ligandStartPos.clear();
        loadListArg(value,currState.ligandStartPos);
    }
    else if(name == "snap_list")
    {
        currState.snapList.clear();
        loadListArg(value,currState.snapList);
    }
    else if(name == "checkpoint")
    {
        currState.checkpointFilename = value;
        currState.checkpoint_in();
    }
    else
    {
      char error[256];
      sprintf(error,"I don't know what to do with the parameter %s (=%s)\n",name.c_str(),value.c_str());
      throw mmpbsa::MMPBSAException(error,mmpbsa::COMMAND_LINE_ERROR);
    }
    return 0;
}

int parseFlag(std::string flag, mmpbsa::MeadInterface& mi)
{
    if(flag.find("=") != flag.npos)
    {
        parseParameter(flag,mi);//Oops that's a flag.
        return 0;
    }

    if(flag == "help" || flag == "h")
      {
	std::cout << helpString() << std::endl;
	return 1;
      }
    else if(flag == "trust_prmtop")
    {
        currState.trustPrmtop = true;
	return 0;
    }

    char error[256];
    sprintf(error,"I don't know what to do with the flag \"%s\"",flag.c_str());
    throw mmpbsa::MMPBSAException(error,mmpbsa::COMMAND_LINE_ERROR);
}

int loadListArg(const std::string& values,std::vector<size_t>& array)
{
    using mmpbsa_utils::StringTokenizer;
    StringTokenizer valTokens(values,",");
    int currValue = 0;
    while(valTokens.hasMoreTokens())
    {
        std::string curr = valTokens.nextToken();
        sscanf(curr.c_str(),"%d",&currValue);
        array.push_back(size_t(currValue));
    }
    return 0;
}

std::string helpString()
{
  return "MMPBSA Calculations\n"
    "Usage: ./mmpbsa prmtop=<prmtop file> mdcrd=<mdcrd file> out=<output file> [optional arguments]\n"
    "\n"
    "radii=<radii file>                   SIZ radii file. If no file is provided, values are used from a \n"
    "          lookup table built into mmpbsa\n"
    "istrength=<strength value>           (default = 0)\n"
    "surf_offset=<surface offset>         (default = 0.92 kcal/mol)\n"
    "surf_tension=<surface tension value> (default = 0.00542 kcal/mol/Ang^2)\n"
    "rec_list=<comma separated list>      (example rec_list=0,13,25)  List of beginning atoms of residues,\n"
    "         where the length is deduced from the parmtop file. Note this is a zero indexed list.\n"
    "lig_list=<comma separated list>      List of beginning atoms of ligand. See also, rec_list.\n"
    "snap_list=<comma separated list>     1-indexed list of snapshots to be include. If this option is not used, all snapshots are calculated."
    "trust_prmtop                         Override the Parmtop sanity check. Use with caution!";
}

MMPBSAState::MMPBSAState()
{
    trustPrmtop = false;
    checkpointFilename = "";
    currentSnap = 0;
    currentMolecule = MMPBSAState::COMPLEX;
    checkpointCounter = 0;
}

bool MMPBSAState::checkpoint_in(const std::string& fileName)
{
    mmpbsa_utils::XMLParser xmlDoc;
    std::string theFileName = fileName;
#ifdef __USE_BOINC__
    int boinc_resolve_filename_s(fileName.c_str(), theFileName);
#endif
    try
    {
        xmlDoc.parse(theFileName);
    }
    catch(mmpbsa::XMLParserException xpe)
    {
        if(this->checkpointCounter == 0)
            return false;
        else
        {
            char error[40+strlen(xpe.what())+fileName.size()];
            xpe.what();
            sprintf(error,"Trouble opening XML file: %s -> %s\n",theFileName.c_str(),xpe.what());
            throw mmpbsa::XMLParserException(error,xpe.getErrType());
        }
    }
    
    std::map<std::string,std::string> checkMap = xmlDoc.getChildren();
    
    bool usedAllParameters = true;
    std::string tag = "";
    for(std::map<std::string,std::string>::iterator it = checkMap.begin();
            it != checkMap.end();it++)
    {
        tag = it->first;
        //chain to load parameters into correct variables.
        if(tag == "current_molecule")
        {
            int currMol = 0;
            sscanf(it->second.c_str(),"%d",&currMol);
            switch(currMol)
            {
                case COMPLEX:
                    this->currentMolecule = COMPLEX;
                    break;
                case LIGAND:
                    this->currentMolecule = LIGAND;
                    break;
                case RECEPTOR:
                    this->currentMolecule = RECEPTOR;
                    break;
                default:
                    this->currentMolecule = END_OF_MOLECULES;
                    usedAllParameters = false;
                    break;
            }
        }//end "current_molecule" case
        else if(tag == "current_snap")
        {
            this->currentSnap = 0;
            sscanf(it->second.c_str(),"%d",&(this->currentSnap));
        }//end "current_snap" case
        else if(tag == "checkpoint_counter")
        {
            this->checkpointCounter = 0;
            sscanf(it->second.c_str(),"%d",&(this->checkpointCounter));
        }
        else
        {
            usedAllParameters = false;
        }
    }
    return usedAllParameters;
}

void MMPBSAState::checkpoint_out()
{
#ifdef __USE_BOINC__
    if(!boinc_time_to_checkpoint())
        return;
#endif
    
    this->checkpointCounter++;
    if(checkpointFilename == "")
        return;

    std::map<std::string,std::string> checkMap;
    char buff[16];
    sprintf(buff,"%d",this->currentMolecule);
    checkMap["current_molecule"] = buff;
    sprintf(buff,"%d",this->currentSnap);
    checkMap["current_snap"] = buff;
    sprintf(buff,"%d",this->checkpointCounter);
    checkMap["checkpoint_counter"] = buff;
    mmpbsa_utils::XMLParser xmlDoc("mmpbsa_state",checkMap);
    xmlDoc.write(this->checkpointFilename);

#ifdef __USE_BOINC__
    boinc_checkpoint_completed();
    if(currState.snapList.size())
        boinc_fraction_done(currState.currentSnap/currState.snapList.size());
#endif
}


void MMPBSAState::flush()
{
    if(outputFile.is_open())
        outputFile.flush();
    if(prmtopFile.is_open())
        prmtopFile.flush();
    if(radiiFile.is_open())
        radiiFile.flush();
    if(trajFile.is_open())
        trajFile.flush();
}

void MMPBSAState::close()
{
    if(outputFile.is_open())
        outputFile.close();
    if(prmtopFile.is_open())
        prmtopFile.close();
    if(radiiFile.is_open())
        radiiFile.close();
    if(trajFile.is_open())
        trajFile.close();
}

