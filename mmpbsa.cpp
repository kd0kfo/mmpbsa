#include "mmpbsa.h"

int main(int argc, char** argv)
{
    try 
    {
        myOutput;
        myOutput.open("/ibis/users_linux/dcoss/working_dir/tut1/polyAT-mmpbsa/myoutput.txt",std::ios::out);
        int returnMe = realDeal(argc, argv);
        myOutput.close();
        return returnMe;
    }    
    catch (MMPBSAException e)
    {
        std::cerr << e.identifier() << ": " << e.what() << std::endl;
        return e.getErrType();
    }
}

int realDeal(int argc, char** argv)
{
    using std::valarray;
    using std::vector;
    using std::slice;
    using std::map;


    MeadInterface mi;

    if(argc == 1)
      {
	parseFlag("help",mi);
	return 0;
      }

    receptorStartPos.push_back(0);//would be better as a command line argument, but I don't know yet what this will look like. To add later.
    ligandStartPos.push_back(1);//array containing the starting positions of ligands and receiptors
    parseArgs(argc,argv,mi);

    //load and check the parmtop file.
    mmpbsa_io::SanderParm * sp = new mmpbsa_io::SanderParm;
    sp->raw_read_amber_parm(::prmtopFile);
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

    for(size_t i = 0;i<receptorStartPos.size();i++)
    {
        size_t currPos = receptorStartPos[i];
        bottom = entireEFun.mol_ranges[2*currPos];
        top = entireEFun.mol_ranges[2*currPos+1];
        valarray<bool> currReceptor(true,top-bottom);
        complexKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorKeepers[slice(bottom,top-bottom,1)] = currReceptor;
        receptorSize += top-bottom;
    }
    for(size_t i = 0;i<ligandStartPos.size();i++)
    {
        size_t currPos = ligandStartPos[i];
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
    if(radiiFile.is_open())
        mmpbsa_io::read_siz_file(radiiFile,radii, residues);

    if(!trajFile.good())
        throw MMPBSAException("Unable to read from trajectory file",BROKEN_TRAJECTORY_FILE);

    using namespace mmpbsa_io;
    get_traj_title(trajFile);//Don't need title, but this ensure we are at the top of the file. If the title is needed later, hook this.
    size_t snapcounter = 0;
    valarray<mmpbsa_t> snapshot(sp->natom*3);
    valarray<mmpbsa_t> complexSnap(complexSize*3);
    valarray<mmpbsa_t> receptorSnap(receptorSize*3);
    valarray<mmpbsa_t> ligandSnap(ligandSize*3);

    while(!trajFile.eof())
    {
        try{
            if(get_next_snap(trajFile, snapshot, sp->natom,true))
                printf("Snapshot #%d has been loaded.\n",++snapcounter);
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

        //if a list of snaps to be run is provided, check to see if this snapshot
        //should be used. Remember: snapcounter is 1-indexed.
        if(::snapList.size())
            if(!mmpbsa_utils::contains(snapList,snapcounter))
            {
                printf("Skipping Snapshot #%d",snapcounter);
                continue;
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

        EMap complexEMap = MeadInterface::full_EMap(complexEFun,complexSnap,fdm,
                *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
        EMap receptorEMap = MeadInterface::full_EMap(receptorEFun,receptorSnap,fdm,
                *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);
        EMap ligandEMap = MeadInterface::full_EMap(ligandEFun,ligandSnap,fdm,
                *pradii,residues,mi.istrength,mi.surf_tension,mi.surf_offset);

        printSnapshot(complexEMap,receptorEMap,ligandEMap,::outputFile);
    }


    delete sp;

    ::outputFile.close();
    ::prmtopFile.close();
    ::radiiFile.close();
    ::trajFile.close();

    return 0;
}

void printSnapshot(const EMap& complexEMap, const EMap& receptorEMap, const EMap& ligandEMap, std::fstream& outFile)
{
    outFile << "COMPLEX" << std::endl;
    outFile << complexEMap  << std::endl;
    outFile << "RECEPTOR" << std::endl;
    outFile << receptorEMap  << std::endl;
    outFile << "LIGAND" << std::endl;
    outFile << ligandEMap  << std::endl << std::endl;
}

void parseArgs(int argc, char** argv, MeadInterface& mi)
{
    using std::string;
    for(int i = 1;i<argc;i++)
    {
        string currArg = argv[i];
        
        if(currArg.substr(0,2) == "--")
            currArg.erase(currArg.begin(),currArg.begin()+1);

        if(currArg.find("=") != string::npos)
            parseParameter(currArg,mi);
        else
            parseFlag(currArg,mi);
    }
}

void parseParameter(std::string arg, MeadInterface& mi)
{
    if(arg.find("=") == string::npos)
        return parseFlag(arg,mi);//oops that's a flag.
    
    if(arg.find("=") != arg.rfind("="))
    {
        char error[256];
        sprintf(error,"Multiple occurance of \"=\" in parameter: %s\n",arg.c_str());
        throw MMPBSAException(error,COMMAND_LINE_ERROR);
    }
    
    using std::string;
    string name = arg.substr(0,arg.find("="));
    string value = arg.substr(arg.find("=")+1);

    if(name == "prmtop" || name == "parmtop")
    {
        ::prmtopFile.open(value.c_str(),std::ios::in);
    }
    else if(name == "coordinates" || name == "traj" || name == "mdcrd")
    {
        ::trajFile.open(value.c_str(),std::ios::in);
    }
    else if(name == "out" || name == "output")
    {
        ::outputFile.open(value.c_str(),std::ios::out);
    }
    else if(name == "radii")
    {
        ::radiiFile.open(value.c_str(),std::ios::in);
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
        ::receptorStartPos.clear();
        loadListArg(value,::receptorStartPos);
    }
    else if(name == "lig_list")
    {
        ::ligandStartPos.clear();
        loadListArg(value,::ligandStartPos);
    }
    else if(name == "snap_list")
    {
        ::snapList.clear();
        loadListArg(value,::snapList);
    }
    else
    {
        fprintf(stderr,"I don't know what to do with the parameter %s (=%s)\n",name.c_str(),value.c_str());
    }
}

void parseFlag(std::string flag, MeadInterface& mi)
{
    if(flag.find("=") != flag.npos)
    {
        parseParameter(flag,mi);//Oops that's a flag.
        return;
    }

    if(flag == "help" || flag == "h")
      {
	std::cout << helpString() << std::endl;
	return;
      }

    fprintf(stderr,"I don't know what to do with the flag \"%s\"\n",flag.c_str());
}

void loadListArg(const std::string& values,std::vector<size_t>& array)
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
    "snap_list=<comma separated list>     1-indexed list of snapshots to be include. If this option is not used, all snapshots are calculated.";
}

