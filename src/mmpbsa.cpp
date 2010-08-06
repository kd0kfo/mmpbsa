#include "mmpbsa.h"


int main(int argc, char** argv)
{
    try 
    {
        ::timeAtPreviousCheckpoint = 0;
        mmpbsa_boinc_init();//must be called before any other BOINC routines. If BOINC is not used, nothing will happen.
        int retval = 0;

        ::processQueue = getQueueFile(argc,argv);
        //If no queue was found, run based on arguments in argv. These should
        //correspond to what is placed in the queue XML file.
        if(processQueue.size() == 0)
        {
            std::map<std::string,std::string> argMap = parseArgs(argc,argv);
            MMPBSAState currState;
            mmpbsa::SanderInterface si;
            retval = parseParameter(argMap,currState,si);
            if(!retval)
            {
                sander_run(currState,si);
                if(!retval && !currState.MDOnly)
                {
                    mmpbsa::MeadInterface mi;
                    retval = parseParameter(argMap,currState,mi);
                    if(!retval)
                        retval = mmpbsa_run(currState,mi);
                }
            }
            else
                retval--;
        }
        else
        {
            int queuePosition= 0;
            for(std::vector<MMPBSAState>::iterator it = processQueue.begin();
                    it != processQueue.end();it++)
            {
                switch(it->currentProcess)
                {
                    case MMPBSAState::SANDER:
                        restart_sander(*it,it->currentSI);
                        if(it->placeInQueue > queuePosition)
                        {
                            queuePosition++;
                            continue;
                        }
                        else if(it->placeInQueue < queuePosition)
                        {
                            it->currentSI.completed = false;
                            it->fractionDone = 0;
                            it->placeInQueue = queuePosition;
                        }
                        retval = sander_run(*it,it->currentSI);
                        break;
                    case MMPBSAState::MMPBSA:
                        restart_mmpbsa(*it);
                        if(it->placeInQueue > queuePosition)
                        {
                            queuePosition++;
                            continue;
                        }
                        else if(it->placeInQueue < queuePosition)
                        {
                            it->fractionDone = 0;
                            it->placeInQueue = queuePosition;
                        }
                        retval = mmpbsa_run(*it,it->currentMI);
                        break;
                }
                if(retval)
                    break;
                queuePosition++;
            }
        }
#ifdef __USE_BOINC__
        if(retval)
            fprintf(stderr,"BOINC Error: %s\n",boincerror(retval));
        boinc_finish(retval);
#endif
        return retval;
    }    
    catch (mmpbsa::MMPBSAException e)
    {
        std::cerr << e.identifier() << ": " << e.what() << std::endl;
#ifdef __USE_BOINC__
        boinc_finish(e.getErrType());
#endif
        return e.getErrType();

    }
}

int mmpbsa_run(MMPBSAState& currState, mmpbsa::MeadInterface& mi)
{
    using std::valarray;
    using std::vector;
    using std::slice;
    using std::map;
    using namespace mmpbsa;

    std::fstream outputFile;
    if(!currState.outputFilename.size())
        currState.outputFilename = "mmpbsa-output.txt";
    outputFile.open(currState.outputFilename.c_str(),std::ios::out | std::ios::app);

    //load and check the parmtop file.
    mmpbsa::SanderParm * sp = new mmpbsa::SanderParm;
    sp->raw_read_amber_parm(currState.prmtopFilename);
    if(!currState.trustPrmtop)
        if(!sp->sanityCheck())
            throw MMPBSAException("Parmtop file is insane.",INVALID_PRMTOP_DATA);

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
    if(currState.radiiFilename.size())
    {
        std::fstream radiiFile(currState.radiiFilename.c_str(),std::ios::in);
        mmpbsa_io::read_siz_file(radiiFile,radii, residues);
        radiiFile.close();
    }

    std::fstream trajFile(currState.trajFilename.c_str(),std::ios::in);
    if(!trajFile.good())
        throw MMPBSAException("Unable to read from trajectory file",BROKEN_TRAJECTORY_FILE);

    using namespace mmpbsa_io;
    get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
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
                mmpbsa_io::skip_next_snap(trajFile,sp->natom,isPeriodic);
                ::updateMMPBSAProgress(currState,1);
                report_boinc_progress();
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
    while(!trajFile.eof())
    {
        try{
            //if a list of snaps to be run is provided, check to see if this snapshot
            //should be used. Remember: snapcounter is 1-indexed.
            snapcounter++;
            if(currState.snapList.size())//check if the current snapshot should be skipped
                if(!mmpbsa_utils::contains(currState.snapList,snapcounter))
                {
                    mmpbsa_io::skip_next_snap(trajFile,sp->natom,isPeriodic);
                    fprintf(stdout,"Skipping Snapshot #%d\n",snapcounter);
                    continue;
                }
        
            if(get_next_snap(trajFile, snapshot, sp->natom,isPeriodic))
                fprintf(stdout,"Running Snapshot #%d\n",snapcounter);
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
            ::updateMMPBSAProgress(currState,1);
            report_boinc_progress();
            continue;//Restarted program at the end of a snapshot. So, move on.
        }
        
        if(currState.currentMolecule == MMPBSAState::COMPLEX)
        {
            EMap complexEMap = MeadInterface::full_EMap(complexEFun,complexSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            outputFile << "COMPLEX" << endl << complexEMap << endl;
            currState.currentMolecule = MMPBSAState::RECEPTOR;
            currState.currentSnap = snapcounter;
            ::updateMMPBSAProgress(currState,0.33333333);
            checkpoint_mmpbsa(currState);
        }

        if(currState.currentMolecule == MMPBSAState::RECEPTOR)
        {
            EMap receptorEMap = MeadInterface::full_EMap(receptorEFun,receptorSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            outputFile << "RECEPTOR" << endl << receptorEMap << endl;
            currState.currentMolecule = MMPBSAState::LIGAND;
            currState.currentSnap = snapcounter;
            ::updateMMPBSAProgress(currState,0.333333);
            checkpoint_mmpbsa(currState);
        }

        if(currState.currentMolecule == MMPBSAState::LIGAND)
        {
            EMap ligandEMap = MeadInterface::full_EMap(ligandEFun,ligandSnap,fdm,
                    *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
            outputFile << "LIGAND" << endl << ligandEMap << endl << endl;
            currState.currentMolecule = MMPBSAState::END_OF_MOLECULES;
            ::updateMMPBSAProgress(currState,0.3333333);
        }
        checkpoint_mmpbsa(currState);
        currState.currentMolecule = MMPBSAState::COMPLEX;//Reset current molecule
    }//end of snapshot loop

    currState.fractionDone = 1.0;
    checkpoint_mmpbsa(currState);

    trajFile.close();
    delete sp;
    return 0;
}

int sander_run(MMPBSAState& currState,mmpbsa::SanderInterface& si)
{
    using namespace mmpbsa;
    if(si.completed || currState.currentProcess == MMPBSAState::MMPBSA)
        return 0;

    int retval = si.start();
    if(retval)
        return retval;

    int status;
    while(true)//si has a fork running. Monitor it.
    {
        if(si.poll(status))
        {
            if(status)
            {
                fprintf(stderr,"Sander(%d) had a problem: 0x%x\n",si.getPID(),status);
                return EXIT_CHILD_FAILED;
            }
            break;
        }
#ifdef __USE_BOINC__
        ::poll_boinc_messages(si);
        if(::boinc_time_to_checkpoint())
            checkpoint_sander(currState,si);
        ::boinc_sleep(SanderInterface::pollPeriod);
#endif
    }
    si.completed = true;
    currState.fractionDone = 1.0;
    checkpoint_sander(currState,si);

    return 0;
}

std::map<std::string,std::string> parseArgs(int argc, char** argv)
{
    using std::string;
    std::map<std::string,std::string> returnMe;

    string name,value;
    for(int i = 1;i<argc;i++)
    {
        string currArg = argv[i];
        int retval = 0;
        if(currArg.substr(0,2) == "--")
            currArg.erase(currArg.begin(),currArg.begin()+2);

        if(currArg.find("=") != string::npos)
        {
            name = currArg.substr(0,currArg.find("="));
            value = currArg.substr(currArg.find("=")+1);
            returnMe[name] = value;
        }
        else
        {
            returnMe[currArg] = "";
        }
    }
    return returnMe;
}

int parseParameter(std::map<std::string,std::string> args, MMPBSAState& currState, mmpbsa::MeadInterface& mi)
{
    for(std::map<std::string,std::string>::const_iterator it = args.begin();it != args.end();it++)
    {
        if (it->second.find("=") != std::string::npos)
        {
            char error[256];
            sprintf(error, "Multiple occurance of \"=\" in parameter: %s=%s\n", it->first.c_str(),it->second.c_str());
            throw mmpbsa::MMPBSAException(error, mmpbsa::COMMAND_LINE_ERROR);
        }

        if (it->first == "istrength")
        {
            mmpbsa_t fValue = 0;
            sscanf(it->second.c_str(), MMPBSA_FORMAT, &fValue);
            mi.istrength = mmpbsa_t(fValue);
        } 
        else if (it->first == "surf_offset")
        {
            mmpbsa_t fValue = 0;
            sscanf(it->second.c_str(), MMPBSA_FORMAT, &fValue);
            mi.surf_offset = mmpbsa_t(fValue);
        } 
        else if (it->first == "surf_tension")
        {
            mmpbsa_t fValue = 0;
            sscanf(it->second.c_str(), MMPBSA_FORMAT, &fValue);
            mi.surf_tension = mmpbsa_t(fValue);
        }
        else if (it->first == "help" || it->first == "h")
        {
            std::cout << helpString() << std::endl;
            return 1;
        }
        else if (it->first == "trust_prmtop")
        {
            currState.trustPrmtop = true;
            return 0;
        }
    }
    return 0;
}

int parseParameter(std::map<std::string,std::string> args, MMPBSAState& currState, mmpbsa::SanderInterface& si)
{
    int returnMe = 0;
    currState.receptorStartPos.push_back(0);//in case these are not set manually by the use. This is the default.
    currState.ligandStartPos.push_back(1);
    for(std::map<std::string,std::string>::const_iterator it = args.begin();it != args.end();it++)
    {
        if (it->second.find("=") != string::npos)
        {
            char error[256];
            sprintf(error, "Multiple occurance of \"=\" in parameter: %s=%s\n", it->first.c_str(),it->second.c_str());
            throw mmpbsa::MMPBSAException(error, mmpbsa::COMMAND_LINE_ERROR);
        }

        if (it->first == "prmtop" || it->first == "parmtop") {
            mmpbsa_io::resolve_filename(it->second, si.prmtopFilename);
            currState.prmtopFilename = si.prmtopFilename;
        } else if (it->first == "inpcrd") {
            mmpbsa_io::resolve_filename(it->second, si.inpcrdFilename);
        } else if (it->first == "mdout") {
            mmpbsa_io::resolve_filename(it->second, si.mdoutFilename);
        } else if (it->first == "in" || it->first == "mdin") {
            mmpbsa_io::resolve_filename(it->second, si.mdinFilename);
        } else if (it->first == "rst" || it->first == "restart" || it->first == "restrt") {
            mmpbsa_io::resolve_filename(it->second, si.restartFilename);
        } else if (it->first == "mdcrd") {
            mmpbsa_io::resolve_filename(it->second, si.mdcrdFilename);
            mmpbsa_io::resolve_filename(it->second, currState.trajFilename);
        } else if (it->first == "checkpoint") {
            mmpbsa_io::resolve_filename(it->second, currState.checkpointFilename);
        } else if (it->first == "radii") {
            mmpbsa_io::resolve_filename(it->second, currState.radiiFilename);
        } else if (it->first == "mmpbsa_out") {
            mmpbsa_io::resolve_filename(it->second, currState.outputFilename);
        } else if (it->first == "rec_list") {
            currState.receptorStartPos.clear();
            loadListArg(it->second, currState.receptorStartPos,1);
        } else if (it->first == "lig_list") {
            currState.ligandStartPos.clear();
            loadListArg(it->second, currState.ligandStartPos,1);
        } else if (it->first == "snap_list" || it->first == "snaplist") {
            currState.snapList.clear();
            loadListArg(it->second, currState.snapList);
        } else if (it->first == "help" || it->first == "h") {
            std::cout << helpString() << std::endl;
            return 1;
        } 
        else if(it->first == "sample_queue")
        {
            if(!it->second.size())
            {
                std::cout << "Parameter \"sample_queue\" requires a file name\"" << std::endl;
                return 2;
            }
            std::string queueFilename;
            mmpbsa_io::resolve_filename(it->second,queueFilename);
            ::sampleQueue(queueFilename);
            return 1;
        }
        else if (it->first == "trust_prmtop") {
            currState.trustPrmtop = true;
            return 0;
        } else if (it->first == "mmpbsa_only") {
            si.completed = true;
            currState.currentProcess = MMPBSAState::MMPBSA;
            return 0;
        } else if (it->first == "md_only") {
            currState.MDOnly = true;
            return 0;
        }
        else if(it->first == "weight")
        {
            sscanf(it->second.c_str(),"%f",&currState.weight);
        }

    }
    return returnMe;
}

int loadListArg(const std::string& values,std::vector<size_t>& array, const size_t& offset)
{
    using mmpbsa_utils::StringTokenizer;
    StringTokenizer valTokens(values,",");
    int currValue = 0;
    while(valTokens.hasMoreTokens())
    {
        std::string curr = valTokens.nextToken();
        sscanf(curr.c_str(),"%d",&currValue);
        array.push_back(size_t(currValue) - offset);
    }
    return 0;
}

std::string helpString()
{
  std::string returnMe = "MMPBSA Calculations\n";
#ifdef __USE_BOINC__
  returnMe += "Compiled with BOINC\n";
#endif

    return returnMe + "Usage: ./mmpbsa prmtop=<prmtop file> mdcrd=<mdcrd file> out=<output file> [optional arguments]\n"
    "\n"
    "radii=<radii file>                   SIZ radii file. If no file is provided, "
    "                                     values are used from a lookup table "
    "                                     built into mmpbsa\n"
    "istrength=<strength value>           (default = 0)\n"
    "surf_offset=<surface offset>         (default = 0.92 kcal/mol)\n"
    "surf_tension=<surface tension value> (default = 0.00542 kcal/mol/Ang^2)\n"
    "rec_list=<comma separated list>      (example rec_list=0,13,25)  List of "
    "                                     beginning atoms of residues,\n"
    "                                     where the length is deduced from the "
    "                                     parmtop file. Note this is a one-indexed list.\n"
    "lig_list=<comma separated list>      List of beginning atoms of ligand. "
    "                                     See also, rec_list.\n"
    "snap_list=<comma separated list>     1-indexed list of snapshots to be "
    "                                     include. If this option is not used, "
    "                                     all snapshots are calculated."
    "trust_prmtop                         Override the Parmtop sanity check. Use with caution!"
    "create_sample=<filename>             Creates a sample queue XML file.";
}

MMPBSAState::MMPBSAState()
{
    trustPrmtop = false;
    checkpointFilename = "";
    currentSnap = 0;
    currentProcess = MMPBSAState::SANDER;
    currentMolecule = MMPBSAState::COMPLEX;
    checkpointCounter = 0;
    this->fractionDone = 0;
    placeInQueue = 0;
    this->weight = 1;
    MDOnly = false;
}

MMPBSAState::MMPBSAState(const MMPBSAState& orig)
{
    receptorStartPos = orig.receptorStartPos;
    ligandStartPos = orig.ligandStartPos;
    snapList = orig.snapList;
    outputFilename = orig.outputFilename;
    prmtopFilename = orig.prmtopFilename;
    trajFilename = orig.trajFilename;
    radiiFilename = orig.radiiFilename;
    checkpointFilename = orig.checkpointFilename;
    checkpointCounter = orig.checkpointCounter;
    currentSnap = orig.currentSnap;
    fractionDone = orig.fractionDone;
    MDOnly = orig.MDOnly;
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
}

MMPBSAState& MMPBSAState::operator=(const MMPBSAState& orig)
{
    if(this == &orig)
        return *this;

    receptorStartPos = orig.receptorStartPos;
    ligandStartPos = orig.ligandStartPos;
    snapList = orig.snapList;
    outputFilename = orig.outputFilename;
    prmtopFilename = orig.prmtopFilename;
    trajFilename = orig.trajFilename;
    radiiFilename = orig.radiiFilename;
    checkpointFilename = orig.checkpointFilename;
    checkpointCounter = orig.checkpointCounter;
    currentSnap = orig.currentSnap;
    fractionDone = orig.fractionDone;
    MDOnly = orig.MDOnly;
    currentSI = orig.currentSI;
    currentMI = orig.currentMI;
    currentMolecule = orig.currentMolecule;
    currentProcess = orig.currentProcess;
    trustPrmtop = orig.trustPrmtop;
    placeInQueue = orig.placeInQueue;
    weight = orig.weight;
    
    return *this;
}

bool restart_sander(MMPBSAState& restartState, mmpbsa::SanderInterface& si)
{
    mmpbsa_utils::XMLParser xmlDoc;
    if(restartState.checkpointFilename.size() == 0)
        return false;

    try{
        xmlDoc.parse(restartState.checkpointFilename);
        if(xmlDoc.mainTag() == "mmpbsa_state")
            restartState.currentProcess = MMPBSAState::MMPBSA;
    }
    catch(mmpbsa::XMLParserException xpe)
    {
        fprintf(stderr,"Did not open %s\n",restartState.checkpointFilename.c_str());
        if(xpe.getErrType() == mmpbsa::FILE_READ_ERROR)
            return false;
        else
            throw xpe;
    }
    std::map<std::string,std::string> checkMap = xmlDoc.getChildren();

    bool usedAllParameters = true;
    std::string tag = "";
    for(std::map<std::string,std::string>::iterator it = checkMap.begin();
            it != checkMap.end();it++)
    {
        tag = it->first;
        //chain to load parameters into correct variables.
        if(tag == "pid")
        {
            int newPID = 0;
            sscanf(it->second.c_str(),"%d",&newPID);
            si.setPID(newPID);
        }
        else if(tag == "cpu_time")
        {
            double new_start_time = 0;
            sscanf(it->second.c_str(),"%e",&new_start_time);
            si.set_start_time(new_start_time);
        }
        else if(tag == "finished")
        {
            int newPID = 0;
            sscanf(it->second.c_str(),"%d",&newPID);
            si.completed = newPID;
        }
        else if(tag == "queue_position")
        {
            sscanf(it->second.c_str(),"%d",&(restartState.placeInQueue));
        }
        else
        {
            usedAllParameters = false;
        }
    }
    return usedAllParameters;
}

bool restart_mmpbsa(MMPBSAState& restartState)
{
    if(restartState.checkpointFilename.size() == 0)
        return false;

    mmpbsa_utils::XMLParser xmlDoc;
    try{
        xmlDoc.parse(restartState.checkpointFilename);
    }
    catch(mmpbsa::XMLParserException xpe)
    {
        if(xpe.getErrType() == mmpbsa::FILE_READ_ERROR)
            return false;
        else
            throw xpe;
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
                case MMPBSAState::COMPLEX:
                    restartState.currentMolecule = MMPBSAState::COMPLEX;
                    break;
                case MMPBSAState::LIGAND:
                    restartState.currentMolecule = MMPBSAState::LIGAND;
                    break;
                case MMPBSAState::RECEPTOR:
                    restartState.currentMolecule = MMPBSAState::RECEPTOR;
                    break;
                default:
                    restartState.currentMolecule = MMPBSAState::END_OF_MOLECULES;
                    usedAllParameters = false;
                    break;
            }
        }//end "current_molecule" case
        else if(tag == "current_snap")
        {
            sscanf(it->second.c_str(),"%d",&(restartState.currentSnap));
        }//end "current_snap" case
        else if(tag == "checkpoint_counter")
        {
            sscanf(it->second.c_str(),"%d",&(restartState.checkpointCounter));
        }
        else if(tag == "queue_position")
        {
            sscanf(it->second.c_str(),"%d",&(restartState.placeInQueue));
        }
        else
        {
            usedAllParameters = false;
        }
    }
    return usedAllParameters;
}

void checkpoint_sander(MMPBSAState& saveState, mmpbsa::SanderInterface& si)
{
    std::fstream sanderProgress("progress.dat",std::ios::in);
    if(sanderProgress.good())
    {
        double completed,remaining;
        char buff[256];
        sanderProgress.getline(buff,256);
        sscanf(buff,"%lf %lf",&completed,&remaining);
        saveState.fractionDone = completed/(completed+remaining);
        sanderProgress.close();
    }
    std::map<std::string,std::string> checkMap;
    char buff[16];
    sprintf(buff,"%d",si.getPID());
    checkMap["pid"] = buff;
    sprintf(buff,"%e",si.start_time());
    checkMap["start_time"] = buff;
    sprintf(buff,"%e",si.cpu_time());
    checkMap["cpu_time"] = buff;
    sprintf(buff,"%e",si.netRuntime());
    checkMap["runtime"] = buff;
    sprintf(buff,"%d",si.completed);
    checkMap["finished"] = buff;
    sprintf(buff,"%d",saveState.placeInQueue);
    checkMap["queue_position"] = buff;
    sprintf(buff,"%lf",saveState.fractionDone);
    checkMap["stage_fraction_done"] = buff;
    mmpbsa_utils::XMLParser xmlDoc("moldyn_state",checkMap);
    checkpoint_out(saveState,xmlDoc);
}

void checkpoint_mmpbsa(MMPBSAState& saveState)
{
    
    std::map<std::string,std::string> checkMap;
    char buff[16];
    sprintf(buff,"%d",saveState.currentMolecule);
    checkMap["current_molecule"] = buff;
    sprintf(buff,"%d",saveState.currentSnap);
    checkMap["current_snap"] = buff;
    sprintf(buff,"%d",saveState.checkpointCounter);
    checkMap["checkpoint_counter"] = buff;
    sprintf(buff,"%d",saveState.placeInQueue);
    checkMap["queue_position"] = buff;
    mmpbsa_utils::XMLParser xmlDoc("mmpbsa_state",checkMap);
    checkpoint_out(saveState,xmlDoc);
}

void checkpoint_out(MMPBSAState& saveState,mmpbsa_utils::XMLParser& xmlDoc)
{
    report_boinc_progress();
    saveState.checkpointCounter++;
    if(saveState.checkpointFilename != "")
        xmlDoc.write(saveState.checkpointFilename);
#ifdef __USE_BOINC__
    boinc_checkpoint_completed();
#endif
}


int mmpbsa_boinc_init()
{
#ifndef __USE_BOINC__
    return 0;
#else
    BOINC_OPTIONS options;
    memset(&options, 0, sizeof(options));
    options.main_program = true;
    options.check_heartbeat = true;
    options.handle_process_control = true;
    fprintf(stderr, "mmpbsa started\n");
    return boinc_init_options(&options);
#endif
}

void send_status_message(mmpbsa::SanderInterface& si, double frac_done, 
        double checkpoint_cpu_time)
{
#ifdef __USE_BOINC__
    double current_cpu_time =  si.start_time() + si.cpu_time();
    boinc_report_app_status(current_cpu_time,checkpoint_cpu_time,frac_done);
#endif
}

void poll_boinc_messages(mmpbsa::SanderInterface& si)
{
#ifdef __USE_BOINC__
    BOINC_STATUS status;
    boinc_get_status(&status);
    if (status.no_heartbeat) {
        si.kill();
        exit(0);
    }
    if (status.quit_request) {
        si.kill();
        exit(0);
    }
    if (status.abort_request) {
        si.kill();
        exit(0);
    }
    if (status.suspended) {
        if (!si.isSuspended()) {
            si.stop();
        }
    } else {
        if (si.isSuspended()) {
            si.resume();
        }
    }
#endif
}


std::vector<MMPBSAState> getQueueFile(int argc,char** argv)
{
    using mmpbsa_utils::XMLParser;
    std::vector<MMPBSAState> returnMe;

    XMLParser queueXML;
    std::string xmlFilename = "";
    std::string arg,name,value;
    for(size_t i = 1;i<argc;i++)
    {
        arg = argv[i];
        if(arg.find("=") == std::string::npos)
            continue;
        name = arg.substr(0,arg.find("="));
        value = arg.substr(arg.find("=")+1);
        if(name == "queue")
        {
            mmpbsa_io::resolve_filename(value,xmlFilename);
            break;
        }
    }

    if(xmlFilename.size() == 0)
        return returnMe;

    try
    {
        queueXML.parse(xmlFilename);
    }
    catch(mmpbsa::XMLParserException xpe)
    {
        fprintf(stderr,"Did not open %s\n",xmlFilename.c_str());
        if(xpe.getErrType() == mmpbsa::FILE_READ_ERROR)
            return returnMe;
        else
            throw xpe;
    }

    xmlNodePtr head = queueXML.getHead();
    if(!head)
        return returnMe;

    if(xmlStrEqual(head->name,(xmlChar*)"grid_queue"))
        head = head->children;

    int queuePosition = 0;
    for(xmlNodePtr sibling = head;sibling;sibling = sibling->next)
    {
        MMPBSAState nodeState;
        mmpbsa::SanderInterface si;
        std::map<std::string,std::string> tags = XMLParser::mapNode(sibling);
        if(xmlStrEqual(sibling->name,(xmlChar*)"mmpbsa"))
        {
            mmpbsa::MeadInterface mi;
            parseParameter(tags,nodeState,si);
            parseParameter(tags,nodeState,mi);
            nodeState.currentMI = mi;
            nodeState.currentProcess = MMPBSAState::MMPBSA;
            nodeState.placeInQueue = queuePosition++;
            returnMe.push_back(nodeState);
        }
        else if(xmlStrEqual(sibling->name,(xmlChar*)"molecular_dynamics"))
        {
            parseParameter(tags,nodeState,si);
            nodeState.currentSI = si;
            nodeState.currentProcess = MMPBSAState::SANDER;
            nodeState.placeInQueue = queuePosition++;
            returnMe.push_back(nodeState);
        }
    }

    return returnMe;
}

void sampleQueue(const std::string& filename)
{
    using std::map;
    using mmpbsa_utils::XMLParser;

    map<std::string,std::string> queueMap;
    queueMap["mdin"] = "sander_input.in";
    queueMap["mdout"] = "sander_output.out";
    queueMap["restart"] = "sander_restart.rst";
    queueMap["inpcrd"] = "sander_input_coordinates.inpcrd";
    queueMap["prmtop"] = "sander_prmtop_file.prmtop";
    queueMap["mdcrd"] = "sander_snapshot_file.mdcrd";
    queueMap["checkpoint"] = "checkpoint_file_name.xml";
    XMLParser sanderXML("molecular_dynamics",queueMap);

    queueMap.clear();

    queueMap["prmtop"] = "sander_prmtop_file.prmtop";
    queueMap["mdcrd"] = "sander_snapshot_file.mdcrd";
    queueMap["radii"] = "DelPhi_radii_file.siz";
    queueMap["mmpbsa_out"] = "mmpbsa-result-output.out";
    queueMap["snap_list"] = "1,3";
    queueMap["checkpoint"] = "checkpoint_file_name.xml";
    XMLParser mmpbsaXML("mmpbsa",queueMap);
    queueMap.clear();

    xmlAddNextSibling(sanderXML.getHead(),mmpbsaXML.getHead());
    XMLParser theDoc("grid_queue",queueMap);
    ::xmlAddChild(theDoc.getHead(),sanderXML.getHead());
    theDoc.write(filename);
}

void updateMMPBSAProgress(MMPBSAState& currState,const double& increment)
{
    if(currState.snapList.size())
    {
        currState.fractionDone += increment/currState.snapList.size();
    }
    else
    {
        currState.fractionDone += increment/5;
    }
}

void report_boinc_progress()
{
#ifdef __USE_BOINC__
    double totalWeight = 0;
    double completed = 0;
    for(size_t i = 0;i< ::processQueue.size();i++)
    {
        totalWeight += processQueue[i].weight;
        completed += processQueue[i].fractionDone*processQueue[i].weight;
    }
    //boinc_fraction_done(completed/totalWeight);
    double cpu_time;
    ::boinc_wu_cpu_time(cpu_time);
    ::boinc_report_app_status(cpu_time,::timeAtPreviousCheckpoint,completed/totalWeight);
#endif
}
